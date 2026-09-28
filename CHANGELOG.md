# Changelog

## Unreleased — 2026-09-28

### Improve end-to-end performance with large BED panels

- Index per-contig BED intervals by prefix-maximum end coordinate and use binary
  search to skip intervals preceding each coverage query. Nested, overlapping,
  and duplicate intervals retain their original counts and report order. Queries
  need not arrive in coordinate order. Unsorted BED input uses the legacy scan.
- Copy the coverage index for post-processing statistics and leave the existing
  UMI clustering, consensus scoring, and coordinate-output rules unchanged.
- Buffer BED JSON output instead of flushing the file after every interval.
  Report values and formatting are unchanged.
- Add coverage regression tests and `benchmarks/compare_outputs.py` to compare
  complete alignment records and JSON statistics between runs.
- All 4,719,377 output alignment records match the pre-optimization output,
  including sequences, qualities, mapping fields and auxiliary tags (ignoring
  order within equal-coordinate groups). All JSON statistics and coverage values
  match after excluding command paths. The output passes quickcheck and indexing;
  26 GoogleTest cases and five end-to-end tests pass.
- On the supplied 6,500,311-record input with the same hg19 reference, WES BED,
  and default two I/O workers, elapsed time decreased from 343.93 to 115.99
  seconds, including JSON/HTML reports (about 2.97x). Peak RSS changed from
  1,997,924 to 2,000,624 KiB. These are individual local WSL runs, not repeated
  controlled trials or a general speed guarantee; storage and caching affect timing.

### Fix output coordinate regression with coordinate-sorted input

**Symptom:** gencore could exit with `ERROR: internal output ordering failure`
even when the input BAM was correctly coordinate-sorted. The partial BAM left
by a failed run must not be treated as a complete result.

**Plain-language explanation:** gencore holds groups of reads in memory while
waiting for their mates and generating consensus reads. To write records in
order, it must know the earliest coordinate still waiting to be processed.
The previous code used the sign of TLEN (template length) to decide which mate
started first, but this assumption did not hold for all alignments. It could
therefore write a later record while an earlier read was still pending.
Writing that earlier read afterwards caused the output coordinates to go backwards.

In the reproduced case, a read started at `chrY:13140999`, its mate started at
`chrY:13141011`, and the first read had TLEN `-131`. The previous code incorrectly
set the group's left boundary to `13141011`, allowing a record at `13141007` to
be written before the pending record at `13140999`.
These are 1-based SAM coordinates; the error uses 0-based coordinates:

```text
ERROR: internal output ordering failure. Found 23:13140998 after 23:13141006
```

Here, `23` is the zero-based index of chrY in this BAM header, not a chromosome
name. The failure was caused by the program's TLEN assumption, rather than
input sorting, UMI processing, or a reference/BED mismatch.

**Changes:**

- Use `min(POS, PNEXT)` as the left boundary in the existing nearby,
  same-contig pairing path, independently of TLEN sign.
- Ensure the finalization boundary covers at least `max(POS, PNEXT)`. For zero
  or short TLEN, wait until input has passed the later mate's coordinate before
  finalizing the group.
- Send records without usable mate coordinates directly to the ordered output
  buffer instead of creating negative-coordinate clusters.
- Retain input-order validation and output-regression errors. The fix does not
  require an additional full-file sort or buffering the entire BAM in memory.

**Compatibility:** No new command-line options are required. Directional UMI
clustering and consensus scoring rules are unchanged. Corrected coordinate
grouping and mate handling can change affected families and consensus results.
Replace the executable and rerun from the original sorted BAM; do not reuse the
partial output from a failed run.

**Validation:**

- Added a minimal reproducer that failed with the previous binary and passes
  after the fix.
- Added regression cases for negative TLEN on the leftmost read, zero TLEN, and
  short TLEN, comparing streaming finalization with end-of-file finalization.
- All 22 GoogleTest cases, five end-to-end tests, and built-in tests passed.
- A complete run using the supplied BAM, hg19 reference, and BED exited with
  status 0. It read 6,500,311 input records and produced 4,719,377 output records.
  The output passed `samtools quickcheck`, was successfully indexed with
  `samtools index` without another sorting step, and was fully readable for counting.
- Rebuilt the fully static Release executable at `bin/gencore`.
