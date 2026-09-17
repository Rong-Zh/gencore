#ifndef GENCORE_H
#define GENCORE_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "options.h"
#include "cluster.h"
#include "stats.h"
#include "bed.h"
#include "htslib/sam.h"
#include <map>
#include <set>
#include <compare>
#include <functional>
#include <optional>
#include "bamutil.h"

using namespace std;

struct BamCoordinate {
    int tid;
    hts_pos_t pos;

    [[nodiscard]] static BamCoordinate from(const bam1_t* b) noexcept {
        return b->core.tid >= 0 && b->core.pos >= 0
            ? BamCoordinate{b->core.tid, b->core.pos}
            : BamCoordinate{-1, -1};
    }

    [[nodiscard]] bool isMapped() const noexcept {
        return tid >= 0;
    }

    std::strong_ordering operator<=>(const BamCoordinate& other) const noexcept {
        if(isMapped() != other.isMapped())
            return isMapped() ? std::strong_ordering::less
                              : std::strong_ordering::greater;
        if(tid != other.tid)
            return tid <=> other.tid;
        return pos <=> other.pos;
    }

    bool operator==(const BamCoordinate&) const noexcept = default;
};

struct bamComp {
    bool operator()(const bam1_t* b1, const bam1_t* b2) const noexcept {
        const auto coordinateOrder = BamCoordinate::from(b1) <=> BamCoordinate::from(b2);
        if(coordinateOrder != 0)
            return coordinateOrder < 0;

        // Coordinate sort does not prescribe the order of ties.  The object
        // address makes distinct records distinct set elements without relying
        // on bam1_t::data or truncating pointers to long.
        return std::less<const bam1_t*>()(b1, b2);
    }
};

class Gencore {
public:
    Gencore(Options *opt);
    ~Gencore();

    void consensus();

private:
	void releaseClusters(map<int, map<int, map<long, Cluster*>>>& clusters);
	void dumpClusters(map<int, map<int, map<long, Cluster*>>>& clusters);
	void addToCluster(bam1_t* b);
	void addToProperCluster(bam1_t* b);
	void addToUnProperCluster(bam1_t* b);
	void createCluster(map<int, map<int, map<long, Cluster*>>>& clusters, int tid, int left, long right);
    void outputPair(Pair* p);
    void finishConsensus(map<int, map<int, map<long, Cluster*>>>& clusters);
    void report();
    void bufferOutput(bam1_t* b);
    void flushReadyOutput(BamCoordinate inputCoordinate);
    void flushOutputBefore(BamCoordinate watermark);
    void outputOutSet();
    void writeBam(bam1_t* b);

private:
    string mInput;
    string mOutput;
    Options *mOptions;
    // chrid:left:right
    map<int, map<int, map<long, Cluster*>>> mProperClusters;
    map<int, map<int, map<long, Cluster*>>> mUnProperClusters;
    bam_hdr_t *mBamHeader;
    samFile* mOutSam;
    hts_tpool* mThreadPool;
    Stats* mPreStats;
    Stats* mPostStats;
    set<bam1_t*, bamComp> mOutSet;
    bool mOutSetCleared;
    bool mProperClustersFinished;
    optional<BamCoordinate> mLastWrittenCoordinate;
    int mClusterTick;
};

#endif
