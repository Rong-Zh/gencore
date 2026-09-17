#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include "umicluster.h"

using Families = std::vector<std::vector<std::string>>;

TEST(UmiClusterTest, MergesSingletonsAndBreaksTiesLexically) {
    EXPECT_EQ(groupUmis({{"AAAT", 1}, {"AAAA", 1}}, 1), (Families{{"AAAA", "AAAT"}}));
}

TEST(UmiClusterTest, KeepsSimilarlyAbundantFamiliesSeparate) {
    EXPECT_EQ(groupUmis({{"AAAA", 10}, {"AAAT", 8}}, 1), (Families{{"AAAA"}, {"AAAT"}}));
}

TEST(UmiClusterTest, UsesDirectionalBoundary) {
    EXPECT_EQ(groupUmis({{"AAAA", 5}, {"AAAT", 3}}, 1), (Families{{"AAAA", "AAAT"}}));
    EXPECT_EQ(groupUmis({{"AAAA", 4}, {"AAAT", 3}}, 1), (Families{{"AAAA"}, {"AAAT"}}));
}

TEST(UmiClusterTest, TraversesErrorChainsBeyondRootDistance) {
    EXPECT_EQ(groupUmis({{"AAAA", 10}, {"AAAT", 3}, {"AATT", 2}}, 1),
              (Families{{"AAAA", "AAAT", "AATT"}}));
}

TEST(UmiClusterTest, DoesNotTreatWeakConnectivityAsOneFamily) {
    EXPECT_EQ(groupUmis({{"AAAA", 10}, {"AAAT", 1}, {"AATT", 8}}, 1),
              (Families{{"AAAA", "AAAT"}, {"AATT"}}));
}

TEST(UmiClusterTest, DoesNotInflateCountsAfterMerging) {
    // Root absorbs AAAT, but 10+4 must not enable absorbing AATA (6).
    EXPECT_EQ(groupUmis({{"AAAA", 10}, {"AAAT", 4}, {"AATA", 6}}, 1),
              (Families{{"AAAA", "AAAT"}, {"AATA"}}));
}

TEST(UmiClusterTest, ZeroThresholdMeansExactMatching) {
    EXPECT_EQ(groupUmis({{"AAAA", 10}, {"AAAT", 1}}, 0), (Families{{"AAAA"}, {"AAAT"}}));
}

TEST(UmiClusterTest, HandlesEmptyInputAndMissingUmis) {
    EXPECT_TRUE(groupUmis({}, 1).empty());
    EXPECT_EQ(groupUmis({{"", 20}}, 1), (Families{{""}}));
    EXPECT_EQ(groupUmis({{"", 20}, {"A", 1}}, 1), (Families{{""}, {"A"}}));
}

TEST(UmiClusterTest, RequiresEqualLengthAndDualUmiStructure) {
    EXPECT_EQ(groupUmis({{"AAAA", 10}, {"AAA", 1}}, 1), (Families{{"AAAA"}, {"AAA"}}));
    EXPECT_EQ(groupUmis({{"AA_A", 10}, {"A_AA", 1}}, 2), (Families{{"AA_A"}, {"A_AA"}}));
    EXPECT_EQ(groupUmis({{"AA_A", 10}, {"AA_T", 1}}, 1), (Families{{"AA_A", "AA_T"}}));
}

TEST(UmiClusterTest, TreatsNAsLiteralNotWildcard) {
    EXPECT_EQ(groupUmis({{"AAAA", 10}, {"AANN", 1}}, 1), (Families{{"AAAA"}, {"AANN"}}));
}

TEST(UmiClusterTest, HandlesLargeCountsWithoutOverflow) {
    const auto maximum = std::numeric_limits<size_t>::max();
    EXPECT_EQ(groupUmis({{"AAAA", maximum}, {"AAAT", maximum}}, 1),
              (Families{{"AAAA"}, {"AAAT"}}));
}

TEST(UmiClusterTest, RejectsNegativeThreshold) {
    EXPECT_THROW(groupUmis({{"AAAA", 1}}, -1), std::invalid_argument);
}
