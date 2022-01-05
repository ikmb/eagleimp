/*
 *    Copyright (C) 2018-2021 by Lars Wienbrandt,
 *    Institute of Clinical Molecular Biology, Kiel University
 *    
 *    This file is part of EagleImp.
 *
 *    EagleImp is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    EagleImp is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with EagleImp. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef TARGET_H_
#define TARGET_H_

#include "VCFData.h"
#include "PBWT.h"

#include "MyMalloc.h"

class Target {
public:
    Target(BooleanVector &phasedTargetMat, BooleanVector &phasedTargetPat,
            vector<float> &phasedDosageMat, vector<float> &phasedDosagePat,
            size_t nt, const string &targetID, const VCFData &vcfdata, size_t NrefhapsAllDiploid, size_t NrefhapsCorrected,
            size_t K, fp_type expectIBDcM, fp_type hist, fp_type pErr, fp_type pLimit,
            bool doPrePhasing, bool doRevPhasing, bool lastIter, bool impMissing, bool skipPhasing, bool useExtArch);

    ~Target() {
        if (bestHapsFlags.getData())
            MyMalloc::free(bestHapsFlags.getData());
        if (monos.getData())
            MyMalloc::free(monos.getData());
        if (splits.getData())
            MyMalloc::free(splits.getData());
    }

    const GenotypeVector& getTargetGTs() const { return targetFull; }

    const BooleanVector& getBestHapsFlags() const { return bestHapsFlags; }

    const BooleanVector& getSplits() const { return splits; }
    size_t getNSplits() const { return splitsites.size(); }
    const vector<size_t>& getSplitSites() const { return splitsites; }

    size_t getIdx() const { return ntarget; }

    // prepPhase() + creates referenceT and inconsT itself (instead of having them provided as below) + phaseExtPBWT()
    void phase();

    // prepare phasing, i.e. prepares required data, e.g. split sites, but NOT the condensed reference
    void prepPhase();

//    // referenceT: condensed haplotype information of the K reference haps at split sites (transposed, i.e. referenceT contains a BooleanVector of length K for each split site)
//    // inconsT: for each haplotype in bestHaps a vector marking each segment (last split site +1 to current split site) with true when it is either inconsistent
//    // to the target genotype or if it is masked out by the IBD check
//    // (transposed, i.e. inconsT contains a BooleanVector of length K for each segment (last split site +1 to current split site))
//    void phaseExtRef(const vector<BooleanVector> &referenceT_, const vector<BooleanVector> &inconsT_);

    // do phasing given a PBWT
    void phaseExtPBWT(PBWT &pbwt);

    int getMissings() const { return totalmissings; }
    int getMonHets() const { return totalmonhets; }

    fp_type getTotalConfidence() const { return totalconf; }
    size_t getNCallSites() const { return ncallsites; }

#ifdef STOPWATCH
    // returns the time of the phasing step resp. the accelerateable part of the phasing step (CPU-only)
    uint64_t getPhaseTime() { return phaseTime; }
    uint64_t getAccPhaseTime() { return accPhaseTime; }
#endif

    static fp_type cMmaxSplit;
//    static size_t MinIBDSplitSitesA;
//    static size_t MinIBDSplitSitesB;
    static size_t HistorySize;
    static size_t HistorySizeFast;
    static size_t BeamWidth;
    static size_t BeamWidthFast;
    static size_t Delta;
    static size_t DeltaFast;
    static fp_type minFixSize;
    static fp_type maxFixSize;
    static fp_type minFixThresh;
    static fp_type maxFixThresh;

private:

    // find the best haplotype references for current target (using sort method)
    void findBestHaps();

    // deprecated!
    // find the best haplotype references for current target (using MinMaxHeap)
// void findBestHapsMinMax();

    // identify singleton sites in references picked by bestHaps
    void findMonos();
    // identify all call locations (het + non-mono) plus split sites and make vectors with condensed information
    void findCallLocs();
    // create condensed transposed reference
    void condensedReferenceT(vector<BooleanVector> &referenceT, vector<int> *count0);
//    // for each reference in bestHaps mark the sites that are inconsistent with the target plus IBD check
//    void findInconsAndCheckIBD(vector<BooleanVector>& inconsT);
    // merge inconsistencies to dest found by comparing genotype g to the haplotype data in src
    void mergeInconsistencies(Genotype g, const BooleanVector::data_type* src, BooleanVector::data_type* dest, size_t size);

    //const VCFData &vcfdata;

    // the phasing destination
    BooleanVector &phasedTargetMat;
    BooleanVector &phasedTargetPat;
    vector<float> &phasedDosageMat;
    vector<float> &phasedDosagePat;
    bool doPrePhasing;
    bool doRevPhasing;
    bool lastIter;
    bool impMissing;
    bool skipPhasing;
    bool useExtArch;

    const vector<BooleanVector> &referenceFull; // full information of reference haps on tgt sites
    const vector<BooleanVector> &referenceFullT; // full information of reference haps on tgt sites (transposed)
    const vector<fp_type> &cMpos;
    size_t NrefhapsAllDiploid; // number of reference haplotypes (if haploid is encoded hom. diploid)
    size_t NrefhapsCorrected; // number of really available reference haplotypes (NrefhapsAllDiploid minus number of haploid samples minus 2 for one target if we are in an iteration > 1)
    size_t MFull; // all sites (not only split sites)
    size_t K; // number of conditioning haplotypes (max. 2*Nref)
    fp_type expectIBDcM;
    fp_type hist; // history factor
    fp_type pErr; // switch error probability
    fp_type pLimit; // probability deviation limit in beam (usually switch error probability or squared sw err prob)

    size_t ntarget; // index of current target
    const GenotypeVector &targetFull; // full genotype information of target

    const VCFData &vcfdata;

    // filled during preparation of phasing:

    // contains indices to references from the reference set
    vector<size_t> bestHaps;
    BooleanVector bestHapsFlags; // length is NrefhapsAD bits, '1' marks a hap as a "best hap"
    // singleton sites (monos) in bestHaps are marked with true
    BooleanVector monos;
    // all split sites (start site of split): the first element is always the first het site in the target, i.e. the first segment is ignored since the phase is clear anyway
    // includes all callLocs + breaks in long homozygous runs
    vector<size_t> splitsites; // the positions of the splits
    BooleanVector splits;      // each site is marked: split (true) or no split (false)
    GenotypeVector target; // condensed genotype information of target at split sites
    // difference between split sites in cM (cMsplitDiffs[1] = cMpos[splits[1]] - cMpos[splits[0]],...) (cMsplitDiffs[0] = cMpos[splits[0]] - cMpos[0]) (cMsplitDiffs[splits.size()] = cMpos.back() - cMpos[splits.back()])
    vector<fp_type> cMsplitDiffs;
    // indices pointing into splits, heterozygous non-mono sites in current target for which phase calls are required
    vector<size_t> callLocs;

    string targetID;

    // average phasing confidence: split in total confidence and call sites
    fp_type totalconf = 0;
    size_t ncallsites = 0;

    int totalmissings;
    int totalmonhets;

#ifdef STOPWATCH
    uint64_t phaseTime = 0;
    uint64_t accPhaseTime = 0;
#endif

};

class ConstraintContainer {
public:
    size_t constraint_idx; // reference!
    fp_type p; // scaled confidence, expected to be between 0.0 and 0.5 (0.0 and 2.0 for dosage confidences (homozygous))
    ConstraintContainer(size_t idx, fp_type p_) : constraint_idx(idx), p(p_) {}
};

#endif /* TARGET_H_ */
