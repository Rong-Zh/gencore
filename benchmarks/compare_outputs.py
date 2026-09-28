"""Compare report statistics and every SAM field, allowing equal-coordinate ties.

Usage: python3 benchmarks/compare_outputs.py old.bam new.bam old.json new.json
Requires samtools on PATH; memory is bounded by one coordinate group plus reports.
"""
import itertools
import json
import subprocess
import sys


def records(process):
    for line in process.stdout:
        fields = line.rstrip('\n').split('\t')
        # Auxiliary tag order is not semantically significant.
        yield fields[:11] + sorted(fields[11:])


def groups(process):
    for key, rows in itertools.groupby(records(process), key=lambda row: (row[2], row[3])):
        yield key, sorted(rows)


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    before = subprocess.Popen(['samtools', 'view', sys.argv[1]], stdout=subprocess.PIPE, text=True)
    after = subprocess.Popen(['samtools', 'view', sys.argv[2]], stdout=subprocess.PIPE, text=True)
    count = 0
    try:
        for left, right in itertools.zip_longest(groups(before), groups(after)):
            if left != right:
                raise RuntimeError(f'BAM difference at {left[0] if left else None} / '
                                   f'{right[0] if right else None}')
            count += len(left[1])
        if before.wait() or after.wait():
            raise RuntimeError('samtools failed while decoding BAM')
    finally:
        for process in (before, after):
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
            process.wait()
    print(f'All {count} alignment records match, including sequence, qualities, coordinates and tags.')
    reports = []
    for path in sys.argv[3:]:
        with open(path) as stream:
            report = json.load(stream)
        report.pop('command', None)  # Different output paths are expected.
        reports.append(report)
    if reports[0] != reports[1]:
        raise RuntimeError('Report statistics differ')
    print('All JSON statistics and coverage values match (excluding command paths).')


if __name__ == '__main__':
    main()
