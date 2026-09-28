# Performance checks

## End-to-end BED coverage and report optimization

Coverage queries now use a per-contig prefix-maximum-end index to skip intervals
that cannot overlap the read. The index handles nested intervals and does not
assume query coordinates increase. It uses O(R) integers for R BED intervals;
each query costs O(log R + K), where K is the remaining scanned candidate span
(deeply nested intervals can still make K large). Unsorted BED input retains
the original scan semantics. BED JSON lines use buffered newlines rather than
one forced file flush per interval.

Using the same 6,500,311-record input, hg19 reference, WES BED and two I/O workers:

| Measurement | Before (6a88412) | After |
| --- | --- | --- |
| Entire program, including JSON/HTML | 343.93 s | 115.99 s |
| User CPU time | 126.78 s | 82.16 s |
| System CPU time | 37.32 s | 15.00 s |
| Peak RSS | 1,997,924 KiB | 2,000,624 KiB |
| Output records | 4,719,377 | 4,719,377 |

These are individual runs on the local WSL/mounted-Windows-drive environment;
the ~2.97x observed speedup is not a portable guarantee or a repeated controlled
benchmark. Reference, alignment processing and both reports are included;
subsequent samtools validation and content comparisons are excluded from timing.
All 4,719,377 alignment records matched, including all auxiliary tags, after
allowing equal-coordinate tie ordering. All JSON statistics and coverage values
also matched, excluding only command paths. The optimized BAM passed quickcheck
and indexing without another sort.

To compare complete output records (including all auxiliary tags) and report
statistics, allowing harmless equal-coordinate ordering differences:

```sh
python3 benchmarks/compare_outputs.py old.bam new.bam old.json new.json
```

This requires samtools. JSON comparison excludes only the recorded command,
whose input/output paths may differ. The comparator buffers one coordinate
group at a time, rather than both entire SAM files, and loads the JSON reports.

## UMI clustering microbenchmark

The distance-one UMI path uses a hash index when a coordinate cluster has more
than 512 distinct nonzero-count UMIs. For each node it enumerates single-base
substitutions (including N), leaving `_` positions intact. Candidates are sorted
back into the original count/lexical order before traversal. Smaller clusters,
distances other than one, and unexpected alphabets retain the scanning path.
No command-line options or directional grouping rules change.

The index holds string views into stable node storage, uses O(U) additional
entries, and is local to one coordinate cluster. Enumeration requires O(U*L)
hash probes; hashing each length-L string adds O(L) work per probe. This removes
the old O(U^2) candidate scan at default distance for large clusters; it is not a
claim that the entire BAM pipeline is linear or faster by the same factor.

Compile and run from the repository root (Linux/WSL):

```sh
g++ -std=c++20 -O3 -Isrc benchmarks/umi_cluster.cpp src/umicluster.cpp -o bin/benchmark_umi
./bin/benchmark_umi
bash benchmarks/check_sanitizers.sh
```

The microbenchmark uses seed 1729, random 14-base dual UMIs, and counts 1..20.
It compares complete ordered families with the frozen implementation from
commit 8a96b90 before reporting timings. Three local WSL runs measured:

| Distinct UMIs | Original time | Optimized time | Speedup |
| --- | --- | --- | --- |
| 256 | 0.29–0.40 ms | 0.26–0.38 ms | 1.05–1.17x |
| 2,048 | 14.0–16.4 ms | 3.75–4.55 ms | 3.56–3.73x |
| 8,192 | 190–219 ms | 16.5–19.4 ms | 9.80–12.91x |

These are synthetic, sparse-UMI microbenchmarks, not whole-BAM speedups; dense
networks, long UMIs, storage and compression can change the result. No real
sample throughput claim is made. Differential tests additionally cover dense
short UMIs, N, dual separators, mixed lengths, ties and distances 0/1/2 across
60 seeded datasets. Existing SAM pipeline tests check sorting and EOF flushing.

The review also reduced vector reallocations in consensus construction, fixed
FR/RR support-count overflow (previously stored in one byte), and replaced
memset over the Stats object with scalar initialization so its vector remains
properly constructed. FR/RR now use HTSlib integer updates, replacing an existing
tag rather than appending a duplicate. Supporting counts above 255/65535 are
covered by regression tests.

Remaining hotspots include pairwise CIGAR containment in consensus selection
and linear duplex-family searches. Their algorithms were not changed in this
pass; representative real BAM profiling is needed to prioritize further work.
