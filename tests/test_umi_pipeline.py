"""Exercise actual EOF/streaming flush paths with coordinate-sorted SAM input."""
import pathlib
import subprocess
import tempfile
import unittest

BINARY = pathlib.Path(__file__).resolve().parents[1] / "bin" / "gencore"


class UmiPipelineTest(unittest.TestCase):
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
