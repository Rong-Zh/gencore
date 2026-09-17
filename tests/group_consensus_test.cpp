#include <gtest/gtest.h>

#include <cstdint>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "bamutil.h"
#include "cluster.h"
#include "options.h"
#include "stats.h"
#include <map>

using BamPtr = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>;

// Required by the reporter objects linked into the test binary.
std::string command;

namespace {

BamPtr makeRead(const std::string& qname, const std::string& umi,
                uint16_t flag, hts_pos_t pos, size_t length) {
    BamPtr read(bam_init1(), &bam_destroy1);
    if(!read)
        return read;

    const uint32_t cigar = bam_cigar_gen(length, BAM_CMATCH);
    const std::string sequence(length, 'A');
    const std::vector<char> qualities(length, 30);
    if(bam_set1(read.get(), qname.size(), qname.c_str(), flag,
                0, pos, 60, 1, &cigar, 0, pos, 0,
                length, sequence.c_str(), qualities.data(), 0) < 0) {
        return BamPtr(nullptr, &bam_destroy1);
    }

    const char miTag[2] = {'M', 'I'};
    const std::string encodedUmi = "UMI_" + umi;
    if(bam_aux_update_str(read.get(), miTag, -1, encodedUmi.c_str()) < 0)
        return BamPtr(nullptr, &bam_destroy1);

    return read;
}

Pair* makePair(Options* options, const std::string& name,
               const std::string& umi, size_t leftLength,
               size_t rightLength) {
    auto left = makeRead(name, umi, BAM_FPAIRED | BAM_FREAD1, 100, leftLength);
    auto right = makeRead(name, umi, BAM_FPAIRED | BAM_FREAD2, 200, rightLength);
    if(!left || !right)
        return nullptr;

    auto pair = std::make_unique<Pair>(options);
    pair->setLeft(left.release());
    pair->setRight(right.release());
    return pair.release();
}

TEST(GroupConsensusTest, UsesOneCanonicalUmiForBothConsensusMates) {
    Options options;
    options.umiPrefix = "UMI";
    Cluster cluster(&options);
    Stats preStats(&options);
    Stats postStats(&options);

    // The shorter left read is selected from family A, while the shorter
    // right read is selected from family B.  Their UMIs differ by one base,
    // which is allowed by the default UMI clustering threshold.
    constexpr const char* umiA = "CCCTTGT_GCACGTT";
    constexpr const char* umiB = "GCCTTGT_GCACGTT";
    Pair* familyA = makePair(&options, "readA:UMI_CCCTTGT_GCACGTT", umiA, 10, 12);
    Pair* familyB = makePair(&options, "readB:UMI_GCCTTGT_GCACGTT", umiB, 12, 10);
    ASSERT_NE(familyA, nullptr);
    ASSERT_NE(familyB, nullptr);

    cluster.addPair(familyA);
    cluster.addPair(familyB);

    auto consensuses = cluster.clusterByUMI(
        options.properReadsUmiDiffThreshold, &preStats, &postStats, false);
    ASSERT_EQ(consensuses.size(), 1U);
    std::unique_ptr<Pair> consensus(consensuses.front());
    ASSERT_NE(consensus, nullptr);
    ASSERT_NE(consensus->mLeft, nullptr);
    ASSERT_NE(consensus->mRight, nullptr);
    EXPECT_EQ(consensus->getUMI(), umiA);
    EXPECT_EQ(BamUtil::getUMI(consensus->mLeft, options.umiPrefix), umiA);
    EXPECT_EQ(BamUtil::getUMI(consensus->mRight, options.umiPrefix), umiA);
}

std::map<std::string, int> consensusFamilies(
    const std::map<std::string, int>& counts, bool reverseOrder = false,
    bool duplex = false) {
    Options options;
    options.umiPrefix = "UMI";
    options.disableDuplex = !duplex;
    Cluster cluster(&options);
    Stats before(&options), after(&options);
    std::vector<std::pair<std::string, int>> ordered(counts.begin(), counts.end());
    if(reverseOrder)
        std::reverse(ordered.begin(), ordered.end());
    int index = 0;
    for(const auto& [umi, count] : ordered) {
        for(int i = 0; i < count; ++i) {
            Pair* pair = makePair(&options, "read" + std::to_string(index++), umi, 10, 10);
            if(!pair) {
                ADD_FAILURE() << "Failed to build synthetic pair";
                return {};
            }
            cluster.addPair(pair);
        }
    }
    std::map<std::string, int> result;
    for(Pair* raw : cluster.clusterByUMI(1, &before, &after, false)) {
        std::unique_ptr<Pair> consensus(raw);
        EXPECT_EQ(BamUtil::getUMI(consensus->mLeft, options.umiPrefix), consensus->getUMI());
        EXPECT_EQ(BamUtil::getUMI(consensus->mRight, options.umiPrefix), consensus->getUMI());
        result.emplace(consensus->getUMI(), consensus->mMergeReads);
    }
    EXPECT_TRUE(cluster.mPairs.empty());
    return result;
}

TEST(GroupConsensusTest, KeepsAbundantNeighborsSeparateDespiteIdenticalSequences) {
    const std::map<std::string, int> counts{{"AAAA", 10}, {"AAAT", 8}};
    EXPECT_EQ(consensusFamilies(counts), counts);
}

TEST(GroupConsensusTest, DirectionalChainUsesOriginalCountsAndCanonicalRoot) {
    const std::map<std::string, int> expected{{"AAAA", 15}};
    const std::map<std::string, int> counts{{"AAAA", 10}, {"AAAT", 3}, {"AATT", 2}};
    EXPECT_EQ(consensusFamilies(counts), expected);
    EXPECT_EQ(consensusFamilies(counts, true), expected);
}

TEST(GroupConsensusTest, DuplexMatchingStillRunsAfterUmiGrouping) {
    // The two reciprocal tags remain distinct SSCS families (distance >1),
    // then the existing duplex stage combines their matching sequences.
    const auto result = consensusFamilies({{"AAAA_CCCC", 2}, {"CCCC_AAAA", 2}}, false, true);
    EXPECT_EQ(result.size(), 1U);
}

} // namespace
