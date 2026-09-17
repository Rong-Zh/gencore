#include <chrono>
#include <iostream>
#include <random>
#include "umicluster.h"
#include "../tests/reference/umi_cluster_quadratic.h"

int main() {
    std::mt19937 random(1729);
    for(size_t size : {32, 256, 2048, 8192}) {
        std::map<std::string, size_t> counts;
        while(counts.size() < size) {
            std::string umi(14, 'A');
            for(char& base : umi)
                base = "ACGT"[random() % 4];
            umi.insert(7, "_");
            counts[umi] = 1 + random() % 20;
        }
        const auto start = std::chrono::steady_clock::now();
        const auto reference = quadraticGroupUmis(counts, 1);
        const auto middle = std::chrono::steady_clock::now();
        const auto optimized = groupUmis(counts, 1);
        const auto end = std::chrono::steady_clock::now();
        if(reference != optimized)
            return 1;
        const double before = std::chrono::duration<double, std::milli>(middle - start).count();
        const double after = std::chrono::duration<double, std::milli>(end - middle).count();
        std::cout << size << " UMIs: reference=" << before << " ms, optimized="
                  << after << " ms, speedup=" << before / after << "x\n";
    }
}
