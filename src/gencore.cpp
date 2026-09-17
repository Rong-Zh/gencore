#include "gencore.h"
#include "bamutil.h"
#include "jsonreporter.h"
#include "htmlreporter.h"

Gencore::Gencore(Options *opt){
    mOptions = opt;
    mBamHeader = NULL;
    mOutSam = NULL;
    mThreadPool = NULL;
    mPreStats = new Stats(opt);
    mPreStats->setPostStats(false);
    mPostStats = new Stats(opt);
    mPostStats->setPostStats(true);
    mOutSetCleared = false;
    mProperClustersFinished = false;
    mClusterTick = 0;
}

Gencore::~Gencore(){
    outputOutSet();
    releaseClusters(mProperClusters);
    releaseClusters(mUnProperClusters);
    if(mBamHeader != NULL) {
        bam_hdr_destroy(mBamHeader);
        mBamHeader = NULL;
    }
    if(mOutSam != NULL) {
        if (sam_close(mOutSam) < 0) {
            cerr << "ERROR: failed to close " << mOutput << endl;
            exit(-1);
        }
        mOutSam = NULL;
    }
    if(mThreadPool != NULL) {
        hts_tpool_destroy(mThreadPool);
        mThreadPool = NULL;
    }
    delete mPreStats;
    delete mPostStats;
}

void Gencore::report() {
    JsonReporter jsonreporter(mOptions);
    jsonreporter.report(mPreStats, mPostStats);
    HtmlReporter htmlreporter(mOptions);
    htmlreporter.report(mPreStats, mPostStats);
}

void Gencore::releaseClusters(map<int, map<int, map<long, Cluster*>>>& clusters) {
    map<int, map<int, map<long, Cluster*>>>::iterator iter1;
    map<int, map<long, Cluster*>>::iterator iter2;
    map<long, Cluster*>::iterator iter3;
    for(iter1 = clusters.begin(); iter1 != clusters.end(); iter1++) {
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++) {
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                delete iter3->second;
            }
        }
    }
}

void Gencore::dumpClusters(map<int, map<int, map<long, Cluster*>>>& clusters) {
    map<int, map<int, map<long, Cluster*>>>::iterator iter1;
    map<int, map<long, Cluster*>>::iterator iter2;
    map<long, Cluster*>::iterator iter3;
    for(iter1 = clusters.begin(); iter1 != clusters.end(); iter1++) {
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++) {
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                 iter3->second->dump();
            }
        }
    }
}

void Gencore::outputOutSet() {
    set<bam1_t*, bamComp>::iterator iter;
    for(iter = mOutSet.begin(); iter!=mOutSet.end(); iter++) {
        writeBam(*iter);
        // delete this bam
        bam_destroy1(*iter);
    }
    mOutSet.clear();
    mOutSetCleared = true;
}

void Gencore::flushOutputBefore(BamCoordinate watermark) {
    while(!mOutSet.empty()) {
        set<bam1_t*, bamComp>::iterator iter = mOutSet.begin();
        bam1_t* b = *iter;
        const BamCoordinate coordinate = BamCoordinate::from(b);

        // Unmapped records sort after all mapped records.  They are retained
        // until the final flush.
        if(!coordinate.isMapped())
            break;
        if(coordinate >= watermark)
            break;

        writeBam(b);
        bam_destroy1(b);
        mOutSet.erase(iter);
    }
}

void Gencore::flushReadyOutput(BamCoordinate inputCoordinate) {
    if(!inputCoordinate.isMapped())
        return;

    // For coordinate-sorted input, an unread alignment cannot start before the
    // current input coordinate.  A consensus record that can start earlier
    // must belong to a cluster already in memory.  Therefore the exclusive
    // output watermark is the earlier of the current input coordinate and the
    // earliest active cluster coordinate.
    BamCoordinate watermark = inputCoordinate;

    if(!mProperClusters.empty()) {
        map<int, map<int, map<long, Cluster*>>>::const_iterator tidIter = mProperClusters.begin();
        if(!tidIter->second.empty()) {
            const BamCoordinate clusterCoordinate{
                tidIter->first,
                tidIter->second.begin()->first
            };
            if(clusterCoordinate < watermark)
                watermark = clusterCoordinate;
        }
    }

    flushOutputBefore(watermark);
}

void Gencore::writeBam(bam1_t* b) {
    //BamUtil::dump(b);
    const BamCoordinate coordinate = BamCoordinate::from(b);
    if(mLastWrittenCoordinate && coordinate < *mLastWrittenCoordinate) {
        cerr << "ERROR: internal output ordering failure. Found "
             << coordinate.tid << ":" << coordinate.pos << " after "
             << mLastWrittenCoordinate->tid << ":" << mLastWrittenCoordinate->pos << endl;
        exit(-1);
    }
    if(sam_write1(mOutSam, mBamHeader, b) <0) {
        error_exit("Writing failed, exiting ...");
    }
    mLastWrittenCoordinate = coordinate;

    mPostStats->addRead(b);
}

void Gencore::bufferOutput(bam1_t* b) {
    pair<set<bam1_t*, bamComp>::iterator,bool> ret = mOutSet.insert(b);
    //cerr << "inserting " << (b)->core.tid << ":" << (b)->core.pos << endl;
    //cerr << "head " << (*mOutSet.begin())->core.tid << ":" << (*mOutSet.begin())->core.pos << endl;
    //cerr << "tail " << (*mOutSet.rbegin())->core.tid << ":" << (*mOutSet.rbegin())->core.pos << endl;
    // pointing to its next
    if(ret.second == false) {
        cerr << "ERROR: attempted to buffer the same BAM record twice" << endl;
        BamUtil::dump(b);
        BamUtil::dump(*ret.first);
        error_exit("Internal output buffer ownership error");
    }
}

void Gencore::outputPair(Pair* p) {
    mPostStats->addMolecule(1, p->mLeft && p->mRight);

    if(mOutSam == NULL || mBamHeader == NULL)
        return ;

    if(p->mLeft) {
        bufferOutput(p->mLeft);
        p->mLeft =  NULL;
    }
    if(p->mRight) {
        bufferOutput(p->mRight);
        // right bam will be put in the mOutSet, so make it NULL to avoid being deleted
        p->mRight =  NULL;
    }
}

void Gencore::consensus(){
    samFile *in;
    in = sam_open(mOptions->input.c_str(), "r");
    if (!in) {
        cerr << "ERROR: failed to open " << mOptions->input << endl;
        exit(-1);
    }

    if(ends_with(mOptions->output, "sam"))
        mOutSam = sam_open(mOptions->output.c_str(), "w");
    else 
        mOutSam = sam_open(mOptions->output.c_str(), "wb");
    if (!mOutSam) {
        cerr << "ERROR: failed to open output " << mOptions->output << endl;
        exit(-1);
    }

    mThreadPool = hts_tpool_init(mOptions->threads);
    if(mThreadPool == NULL) {
        error_exit("Failed to create HTSlib thread pool");
    }
    htsThreadPool ioThreadPool{
        .pool = mThreadPool,
        .qsize = mOptions->threads * 2
    };
    if(hts_set_thread_pool(in, &ioThreadPool) < 0 ||
       hts_set_thread_pool(mOutSam, &ioThreadPool) < 0) {
        error_exit("Failed to configure HTSlib thread pool");
    }

    mBamHeader = sam_hdr_read(in);
    mOptions->bamHeader = mBamHeader;
    mPreStats->makeGenomeDepthBuf();
    mPreStats->makeBedStats();
    mPostStats->makeGenomeDepthBuf();
    mPostStats->makeBedStats(mPreStats->mBedStats);

    if (mBamHeader == NULL || mBamHeader->n_targets == 0) {
        cerr << "ERROR: this SAM file has no header " << mInput << endl;
        exit(-1);
    }
    BamUtil::dumpHeader(mBamHeader);

    // The input is verified as coordinate sorted while it is read, and the
    // output path below preserves that order.  Keep the header consistent even
    // when an otherwise valid input header omitted or misstated SO.
    if(sam_hdr_update_hd(mBamHeader, "SO", "coordinate") < 0) {
        cerr << "failed to update output sort order in header" << endl;
        exit(-1);
    }

    if (sam_hdr_write(mOutSam, mBamHeader) < 0) {
        cerr << "failed to write header" << endl;
        exit(-1);
    }

    bam1_t *b = NULL;
    b = bam_init1();
    int r;
    int count = 0;
    bool hasPE = false;
    bool isFirst = true;
    optional<BamCoordinate> lastInputCoordinate;
    while ((r = sam_read1(in, mBamHeader, b)) >= 0) {
        // for the first read, check UMI prefix automatically
        if(isFirst) {
            if(mOptions->umiPrefix == "auto") {
                string umi = BamUtil::getQName(b);
                if(umi.find("umi_") != string::npos) 
                    mOptions->umiPrefix = "umi";
                else if(umi.find("UMI_") != string::npos)
                    mOptions->umiPrefix = "UMI";
                else
                    mOptions->umiPrefix = "";

                if(!mOptions->umiPrefix.empty())
                    cerr << endl << "Detected UMI prefix: " << mOptions->umiPrefix << endl << endl;
            }
            isFirst = false;
        }
        mPreStats->addRead(b);
        count++;
        if(count < 1000) {
            if(b->core.mtid >= 0)
                hasPE = true;
        }
        if(count == 1000 && hasPE == false) {
            cerr << "WARNING: seems that the input data is single-end, gencore will not make consensus read and remove duplication for SE data since grouping by coordination will be inaccurate." << endl << endl;
        }

        // Check coordinate order.  In a coordinate-sorted BAM, unmapped
        // records form the final section, so a mapped record must never occur
        // after one of them.
        const BamCoordinate inputCoordinate = BamCoordinate::from(b);
        const bool isUnmapped = !inputCoordinate.isMapped();
        if(lastInputCoordinate && inputCoordinate < *lastInputCoordinate) {
            cerr << "ERROR: the input is unsorted. Found "
                 << inputCoordinate.tid << ":" << inputCoordinate.pos << " after "
                 << lastInputCoordinate->tid << ":" << lastInputCoordinate->pos << endl;
            cerr << "Please sort the input first." << endl << endl;
            BamUtil::dump(b);
            exit(-1);
        }
        // for testing, we only process to some contig
        if(mOptions->maxContig>0 && b->core.tid>=mOptions->maxContig){
            break;
        }
        // if debug flag is enabled, show which contig we are start to process
        if(mOptions->debug && inputCoordinate.isMapped() &&
           (!lastInputCoordinate || inputCoordinate.tid > lastInputCoordinate->tid)) {
            cerr << "Starting contig " << b->core.tid << endl;
        }
        lastInputCoordinate = inputCoordinate;
        // unmapped reads, we just write it and continue
        if(isUnmapped) {
            // we arrived the end of bam file with unmapped reads, go clear the output set first
            if(!mOutSetCleared) {
                if(!mProperClustersFinished) {
                    mProperClustersFinished = true;
                    finishConsensus(mProperClusters);
                }
                outputOutSet();
            }
            //writeBam(b);
            continue;
        }

        // for secondary alignments, we just skip it
        if(!BamUtil::isPrimary(b)) {
            flushReadyOutput(inputCoordinate);
            continue;
        }
        addToCluster(b);
        flushReadyOutput(inputCoordinate);
        b = bam_init1();
    }

    if(r < -1) {
        error_exit("Failed while reading input BAM/SAM");
    }

    if(!mProperClustersFinished) {
        mProperClustersFinished = true;
        finishConsensus(mProperClusters);
    }
    outputOutSet();
    
    //finishConsensus(mUnProperClusters);

    bam_destroy1(b);
    sam_close(in);

    cerr << "----Before gencore processing:" << endl;
    mPreStats->print();

    cerr << endl << "----After gencore processing:" << endl;
    mPostStats->print();

    report();
}

void Gencore::addToProperCluster(bam1_t* b) {
    int tid = b->core.tid;
    int left = b->core.pos;
    long right;

    if(b->core.mtid == b->core.tid && abs(b->core.mpos - b->core.pos) < 100000) { // process pair synchronously when they are on same contig without huge gap
        if(b->core.isize < 0) {
            left = b->core.mpos;
        }
        right = left + abs(b->core.isize) -  1;
    } else { // cross contig, we only process this read, but dont process its mate
        // no mate or mate is not mapped, we cannot remove duplication or make consensus read, so just write it
        if(b->core.mtid < 0) {
            bufferOutput(b);
            return;
        } else { // cross contig pair mapping
            right = -1L * (long)mBamHeader->target_len[b->core.tid] * (long)(b->core.mtid+1) + (long)b->core.mpos;
        }
    }

    createCluster(mProperClusters, tid, left, right);
    mProperClusters[tid][left][right]->addRead(b);


    mClusterTick++;
    if(mClusterTick % 10000 != 0)
        return;

    // make consensus merge
    map<int, map<int, map<long, Cluster*>>>::iterator iter1;
    map<int, map<long, Cluster*>>::iterator iter2;
    map<long, Cluster*>::iterator iter3;
    bool needBreak = false;
    for(iter1 = mProperClusters.begin(); iter1 != mProperClusters.end();) {
        if(iter1->first > tid || needBreak) {
            break;
        }
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ) {
            if(iter1->first == tid && iter2->first >= b->core.pos) {
                needBreak = true;
                break;
            }
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ) {
                // only deal with the clusters with right < processing pos
                if(iter1->first == tid && iter3->first >= b->core.pos) {
                    break;
                }
                vector<Pair*> csPairs = iter3->second->clusterByUMI(mOptions->properReadsUmiDiffThreshold, mPreStats, mPostStats, iter3->first < 0);
                for(size_t i=0; i<csPairs.size(); i++) {
                    //csPairs[i]->dump();
                    outputPair(csPairs[i]);
                    delete csPairs[i];
                }
                // this tid:left:right is done
                delete iter3->second;
                iter3 = iter2->second.erase(iter3);
            }
            // this tid:left is done
            if(iter2->second.size() == 0) {
                iter2 = iter1->second.erase(iter2);
            } else {
                iter2++;
            }
        }
        // this tid is done
        if(iter1->second.size() == 0) {
            iter1 = mProperClusters.erase(iter1);
        } else {
            iter1++;
        }
    }
}

void Gencore::finishConsensus(map<int, map<int, map<long, Cluster*>>>& clusters) {
    // make consensus merge
    map<int, map<int, map<long, Cluster*>>>::iterator iter1;
    map<int, map<long, Cluster*>>::iterator iter2;
    map<long, Cluster*>::iterator iter3;
    for(iter1 = clusters.begin(); iter1 != clusters.end();) {
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ) {
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ) {
                // for unmapped reads, we just store them
                if(iter1->first < 0 || iter2->first < 0 ) {
                    map<string, Pair*>::iterator iterOfPairs;
                    for(iterOfPairs = iter3->second->mPairs.begin(); iterOfPairs!=iter3->second->mPairs.end(); iterOfPairs++) {
                        //csPairs[i]->dump();
                        outputPair(iterOfPairs->second);
                        delete iterOfPairs->second;
                    }
                } else {
                    const int umiThreshold = &clusters == &mProperClusters
                        ? mOptions->properReadsUmiDiffThreshold
                        : mOptions->unproperReadsUmiDiffThreshold;
                    vector<Pair*> csPairs = iter3->second->clusterByUMI(umiThreshold, mPreStats, mPostStats, iter3->first < 0);
                    for(size_t i=0; i<csPairs.size(); i++) {
                        //csPairs[i]->dump();
                        outputPair(csPairs[i]);
                        delete csPairs[i];
                    }
                }
                // this tid:left:right is done
                delete iter3->second;
                iter3 = iter2->second.erase(iter3);
            }
            // this tid:left is done
            if(iter2->second.size() == 0) {
                iter2 = iter1->second.erase(iter2);
            } else {
                iter2++;
            }
        }
        // this tid is done
        if(iter1->second.size() == 0) {
            iter1 = clusters.erase(iter1);
        } else {
            iter1++;
        }
    }
}

void Gencore::addToUnProperCluster(bam1_t* b) {
    int tid = b->core.tid;
    int left = b->core.pos;
    long right = b->core.mpos;
    if(b->core.mtid < b->core.tid) {
        tid = b->core.mtid;
        left = b->core.mpos;
        right = b->core.pos;
    }
    createCluster(mUnProperClusters, tid, left, right);
    mUnProperClusters[tid][left][right]->addRead(b);
}

void Gencore::createCluster(map<int, map<int, map<long, Cluster*>>>& clusters, int tid, int left, long right) {
    map<int, map<int, map<long, Cluster*>>>::iterator iter1 = clusters.find(tid);

    if(iter1 == clusters.end()) {
        clusters[tid] = map<int, map<long, Cluster*>>();
        clusters[tid][left] = map<long, Cluster*>();
        clusters[tid][left][right] = new Cluster(mOptions);
    } else {
        map<int, map<long, Cluster*>>::iterator iter2  =iter1->second.find(left);
        if(iter2 == iter1->second.end()) {
            clusters[tid][left] = map<long, Cluster*>();
            clusters[tid][left][right] = new Cluster(mOptions);
        } else {
            map<long, Cluster*>::iterator iter3 = iter2->second.find(right);
            if(iter3 == iter2->second.end())
                clusters[tid][left][right] = new Cluster(mOptions);
        }
    }
}

void Gencore::addToCluster(bam1_t* b) {
    // unproperly mapped
    if(b->core.tid < 0) {
        // actually this will never happen since it would be written directly if it's unmapped
        addToUnProperCluster(b);
    } else {
        addToProperCluster(b);
    }
}
