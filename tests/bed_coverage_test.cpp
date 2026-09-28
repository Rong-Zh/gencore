#include <gtest/gtest.h>
#include <algorithm>
#include <fstream>
#include <memory>
#include <random>
#include <unistd.h>
#include "bed.h"

namespace {
class BedCoverageTest : public ::testing::Test {
protected:
    Options options;
    std::string path;
    void SetUp() override {
        char pattern[] = "/tmp/gencore-bed-XXXXXX";
        const int descriptor = mkstemp(pattern);
        ASSERT_GE(descriptor, 0);
        close(descriptor);
        path = pattern;
        options.bedFile = path;
        const std::string header = "@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n";
        options.bamHeader = sam_hdr_parse(header.size(), header.c_str());
        ASSERT_NE(options.bamHeader, nullptr);
    }
    void TearDown() override {
        sam_hdr_destroy(options.bamHeader);
        if(!path.empty())
            unlink(path.c_str());
    }
    void compare(const std::vector<std::pair<int, int>>& intervals) {
        {
            std::ofstream file(path);
            for(const auto& [start, end] : intervals)
                file << "chr1\t" << start << '\t' << end << "\ttarget\n";
        }
        Bed bed(&options);
        bed.loadFromFile();
        Bed copied(&options);
        copied.copyFrom(&bed);
        std::vector<long> expected(intervals.size(), 0);
        std::mt19937 random(37);
        for(int iteration = 0; iteration < 5000; ++iteration) {
            const int start = random() % 11000;
            const int length = random() % 300;
            const int end = start + length;
            // Queries deliberately arrive out of order, as they can for
            // overlapping CIGAR blocks or independently emitted consensus.
            bed.statDepth(0, start, length);
            copied.statDepth(0, start, length);
            for(size_t p = 0; p < intervals.size(); ++p) {
                if(intervals[p].second < start)
                    continue;
                if(intervals[p].first > end)
                    break;
                expected[p] += std::min(intervals[p].second, end) - std::max(intervals[p].first, start);
            }
        }
        bed.statDepth(-1, 0, 100);
        bed.statDepth(1, 0, 100); // Empty contig.
        for(size_t p = 0; p < expected.size(); ++p) {
            EXPECT_EQ(bed.mContigRegions[0][p].mCount, expected[p]) << p;
            EXPECT_EQ(copied.mContigRegions[0][p].mCount, expected[p]) << p;
        }
    }
};

TEST_F(BedCoverageTest, MatchesLegacyForOverlappingNestedAndDuplicateIntervals) {
    std::vector<std::pair<int, int>> intervals{{0, 10000}, {0, 5}, {50, 50}, {50, 150}, {50, 150}};
    for(int start = 100; start < 10000; start += 13)
        intervals.emplace_back(start, start + 30);
    compare(intervals);
}

TEST_F(BedCoverageTest, MatchesLegacyForDisjointIntervalsAndGaps) {
    std::vector<std::pair<int, int>> intervals;
    for(int start = 0; start < 10000; start += 50)
        intervals.emplace_back(start, start + 10);
    compare(intervals);
}

TEST_F(BedCoverageTest, UnsortedInputRetainsLegacyBehavior) {
    compare({{1000, 1100}, {100, 150}, {50, 10000}, {2000, 2100}});
}

TEST_F(BedCoverageTest, EmptyBed) {
    compare({});
}
} // namespace
