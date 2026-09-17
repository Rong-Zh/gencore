#include <gtest/gtest.h>
#include <memory>
#include "bamutil.h"

namespace {
std::string fromTag(const std::string& value, const std::string& prefix) {
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read(bam_init1(), bam_destroy1);
    if(!read || bam_set1(read.get(), 4, "read", BAM_FUNMAP, -1, -1, 0,
                         0, nullptr, -1, -1, 0, 0, nullptr, nullptr, 0) < 0)
        return "FAILED";
    if(bam_aux_update_str(read.get(), "MI", -1, value.c_str()) < 0)
        return "FAILED";
    return BamUtil::getUMI(read.get(), prefix);
}

TEST(UmiParseTest, ParsesActualFastpNameAndDocumentedForms) {
    EXPECT_EQ(BamUtil::getUMI("E250200459L1C036R01600126739:UMI_CCCTTGT_GCACGTT", "UMI"),
              "CCCTTGT_GCACGTT");
    EXPECT_EQ(BamUtil::getUMI("read:UMI_ACGN", "UMI"), "ACGN");
    EXPECT_EQ(BamUtil::getUMI("read:ACGT_TGCA", ""), "ACGT_TGCA");
    EXPECT_EQ(BamUtil::getUMI("read/1:UMI_ACGT_TGCA /1", "UMI"), "ACGT_TGCA");
    EXPECT_EQ(BamUtil::getUMI("read:UMI_ACGT_TGCA/2", "UMI"), "ACGT_TGCA");
    EXPECT_EQ(BamUtil::getUMI("read:umi_ACGT_TGCA", "umi"), "ACGT_TGCA");
    EXPECT_EQ(BamUtil::getUMI("read:TAG_ACGT_TGCA", "TAG"), "ACGT_TGCA");
}

TEST(UmiParseTest, RecognizesTagFormatBeforeRawPayload) {
    EXPECT_EQ(fromTag("UMI_ACGT_TGCA", "UMI"), "ACGT_TGCA");
    EXPECT_EQ(fromTag("read:UMI_ACGT_TGCA", "UMI"), "ACGT_TGCA");
    EXPECT_EQ(fromTag("ACGT_TGCA", "UMI"), "ACGT_TGCA");
    EXPECT_EQ(fromTag("AT_CG", "AT"), "CG");
    EXPECT_EQ(fromTag("AT_CG", "UMI"), "AT_CG");
}

TEST(UmiParseTest, RejectsEmptySegmentsAndPartialPayloads) {
    for(const std::string payload : {"_", "_ACGT", "ACGT_", "AC__GT", "A_C_G", "ACGTX", ""}) {
        EXPECT_EQ(fromTag(payload, "UMI"), "") << payload;
        EXPECT_EQ(fromTag("UMI_" + payload, "UMI"), "") << payload;
        EXPECT_EQ(BamUtil::getUMI("read:UMI_" + payload, "UMI"), "") << payload;
    }
}

TEST(UmiParseTest, DoesNotFindPrefixInsideAnotherField) {
    EXPECT_EQ(BamUtil::getUMI("read:NOTUMI_ACGT", "UMI"), "");
    EXPECT_EQ(BamUtil::getUMI("read:UMI_ACGT:other", "UMI"), "");
    EXPECT_EQ(BamUtil::getUMI("read comment:UMI_ACGT", "UMI"), "");
}
}
