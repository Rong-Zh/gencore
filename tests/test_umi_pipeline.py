"""Exercise actual EOF/streaming flush paths with coordinate-sorted SAM input."""
import pathlib
import subprocess
import tempfile
import unittest

BINARY = pathlib.Path(__file__).resolve().parents[1] / "bin" / "gencore"


class UmiPipelineTest(unittest.TestCase):
    def boundary_case(self, first_tlen, second_tlen, mate_pos, padding=True):
        """Put a flush between two mates while another cluster can be emitted."""
        with tempfile.TemporaryDirectory(prefix="gencore-boundary-") as directory:
            work = pathlib.Path(directory)
            reference = work / "ref.fa"
            reference.write_text(">chr1\n" + "A" * 1000 + "\n>chr2\n" + "A" * 1000 + "\n")
            def row(name, flag, pos, mate_ref, mate, tlen):
                return (pos, f"{name}\t{flag}\tchr1\t{pos}\t60\t10M\t{mate_ref}"
                             f"\t{mate}\t{tlen}\tAAAAAAAAAA\tIIIIIIIIII\n")
            records = [row("target:UMI_AAAA", 83, 101, "=", mate_pos, first_tlen),
                       row("target:UMI_AAAA", 163, mate_pos, "=", 101, second_tlen)]
            if padding:
                # First read + 9,998 records + trigger = 10,000 input records.
                # The cross-contig family at 109 can finish before the mate.
                records.extend(row(f"pad{i}:UMI_CCCC", 65, 109, "chr2", 500, 0)
                               for i in range(9998))
                records.append(row("trigger:UMI_GGGG", 65, 111, "chr2", 600, 0))
            source = work / "input.sam"
            source.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
                              "@SQ\tSN:chr2\tLN:1000\n"
                              + "".join(record for _, record in sorted(records)))
            output = work / "output.sam"
            result = subprocess.run([str(BINARY), "-i", str(source), "-o", str(output),
                                     "-r", str(reference), "--umi_prefix", "UMI",
                                     "--no_duplex", "-@", "1"],
                                    cwd=work, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = [line.split("\t") for line in output.read_text().splitlines()
                    if not line.startswith("@")]
            positions = [int(row[3]) for row in rows]
            self.assertEqual(positions, sorted(positions))
            target = [row for row in rows if row[0] == "target:UMI_AAAA"]
            self.assertEqual([int(row[3]) for row in target], [101, mate_pos])
            return target

    def test_negative_tlen_on_leftmost_read_does_not_advance_watermark(self):
        self.assertEqual(self.boundary_case(-131, 131, 113),
                         self.boundary_case(-131, 131, 113, padding=False))

    def test_zero_tlen_waits_for_later_mate(self):
        self.assertEqual(self.boundary_case(0, 0, 201),
                         self.boundary_case(0, 0, 201, padding=False))

    def test_short_tlen_waits_for_later_mate(self):
        self.assertEqual(self.boundary_case(5, -5, 201),
                         self.boundary_case(5, -5, 201, padding=False))

    def run_case(self, padding, distance):
        with tempfile.TemporaryDirectory(prefix="gencore-umi-") as directory:
            work = pathlib.Path(directory)
            reference = work / "ref.fa"
            reference.write_text(">chr1\n" + "A" * 200000 + "\n")
            records = []

            def pair(name, umi, left, right):
                length = right - left + 10
                for flag, pos, mate, tlen in ((99, left, right, length),
                                               (147, right, left, -length)):
                    records.append((pos, f"{name}:UMI_{umi}\t{flag}\tchr1\t{pos}\t60\t10M"
                                    f"\t=\t{mate}\t{tlen}\tAAAAAAAAAA\tIIIIIIIIII\n"))

            pair("first", "AAAA", 101, 201)
            pair("second", "AAAT", 101, 201)
            # >10,000 records force the first cluster to flush during reading.
            for i in range(padding):
                pair(f"padding{i}", "CCCC", 1001 + i * 30, 1016 + i * 30)
            source = work / "input.sam"
            source.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:200000\n"
                              + "".join(record for _, record in sorted(records)))
            output = work / "output.sam"
            result = subprocess.run(
                [str(BINARY), "-i", str(source), "-o", str(output), "-r", str(reference),
                 "--no_duplex", "-d", str(distance), "-@", "1"],
                cwd=work, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = [line.split("\t") for line in output.read_text().splitlines()
                    if not line.startswith("@")]
            positions = [int(row[3]) for row in rows]
            self.assertEqual(positions, sorted(positions))
            target = [row for row in rows if int(row[3]) < 1000]
            self.assertEqual(len(target), 2 if distance == 1 else 4)
            if distance == 1:
                self.assertTrue(all("MI:Z:AAAA" in row[11:] for row in target))
            return [(row[3], row[9], tuple(tag for tag in row[11:] if tag.startswith("MI:")))
                    for row in target]

    def test_directional_is_identical_at_eof_and_streaming_flush(self):
        self.assertEqual(self.run_case(0, 1), self.run_case(5000, 1))

    def test_zero_distance_keeps_distinct_umis_at_eof(self):
        self.run_case(0, 0)


if __name__ == "__main__":
    unittest.main()
