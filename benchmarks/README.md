# Performance checks

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
