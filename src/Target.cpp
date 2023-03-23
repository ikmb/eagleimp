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

#include <algorithm>
#include <iomanip>
#include <sys/resource.h>
#include <sys/time.h>

#include "MinMaxHeap.h"
#include "BooleanHistory.h"
#include "Datatypes.h"
#include "Beam.h"
#include "utils.h"
#include "StatusFile.h"

#include "Target.h"

using namespace std;

bool operator<(const ConstraintContainer& a, const ConstraintContainer& b) { return a.p < b.p; }

// all defaults. overwritten by main from user options.
/* static */ fp_type Target::cMmaxSplit = 0.5;
///* static */ size_t Target::MinIBDSplitSitesA = 40;
///* static */ size_t Target::MinIBDSplitSitesB = 20;
/* static */ size_t Target::HistorySize = 100;    // must be greater than Delta+1
/* static */ size_t Target::HistorySizeFast = 30; // must be greater than Delta+2
/* static */ size_t Target::BeamWidth = 50; //200;
/* static */ size_t Target::BeamWidthFast = 30; //120;
/* static */ size_t Target::Delta = 20;    // corresponds to callLength=Delta+1 in Eagle. Mustn't be greater than 31!!
/* static */ size_t Target::DeltaFast = 10; // corresponds to callLength=Delta+1 in Eagle. Mustn't be greater than 31!!
/* static */ fp_type Target::minFixSize = 0.5; // at least 50% have to be constrained after the first run
/* static */ fp_type Target::maxFixSize = 0.9; // at most 90% should be constrained (if not 100% confident)
/* static */ fp_type Target::minFixThresh = 0.99; // corresponds to 99% confidence
/* static */ fp_type Target::maxFixThresh = 1.0; // corresponds to 100% confidence

Target::Target(BooleanVector &phasedTargetMat_, BooleanVector &phasedTargetPat_,
        vector<float> &phasedDosageMat_, vector<float> &phasedDosagePat_,
        size_t ntarget_, const string &targetID_, const VCFData &vcfdata_, size_t NrefhapsAllDiploid_, size_t NrefhapsCorrected_,
        size_t K_, fp_type expectIBDcM_, fp_type hist_, fp_type pErr_, fp_type pLimit_,
        bool doPrePhasing_, bool doRevPhasing_, bool lastIter_, bool impMissing_, bool skipPhasing_, bool useExtArch_)
    : phasedTargetMat(phasedTargetMat_),
      phasedTargetPat(phasedTargetPat_),
      phasedDosageMat(phasedDosageMat_),
      phasedDosagePat(phasedDosagePat_),
      doPrePhasing(doPrePhasing_),
      doRevPhasing(doRevPhasing_),
      lastIter(lastIter_),
      impMissing(impMissing_),
      skipPhasing(skipPhasing_),
      useExtArch(useExtArch_),
      referenceFull(vcfdata_.getReference_const()),
      referenceFullT(vcfdata_.getReferenceT_const()),
      cMpos(vcfdata_.getCMPos()),
      NrefhapsAllDiploid(NrefhapsAllDiploid_),
      NrefhapsCorrected(NrefhapsCorrected_),
      MFull(vcfdata_.getNSNPs()),
      K(K_),
      expectIBDcM(expectIBDcM_),
      hist(hist_),
      pErr(pErr_),
      pLimit(pLimit_),
      ntarget(ntarget_),
      targetFull(vcfdata_.getTargets()[ntarget_]),
      vcfdata(vcfdata_),
      targetID(targetID_),
      totalconf(0.0),
      ncallsites(0),
      totalmissings(0),
      totalmonhets(0)
{}


void Target::prepPhase() {

#ifdef DEBUG_TARGET
    cout << "\n" << ntarget << ": " << targetID << "\n" << MFull << " sites." << endl;
#endif

    // find best fitting references
#ifdef STOPWATCH
    Stopwatch swbhs("BestHapsSort");
#endif
    findBestHaps();
#ifdef STOPWATCH
    swbhs.stop();
#endif

#ifdef DEBUG_TARGET
#ifdef DEBUG_PRINT_BESTHAPS
    if (ntarget == 0) {
        cout << "Best haps:" << endl;
//        for (size_t i = 0; i < K; i++)
        for (size_t i = 0; i < 10; i++)
            cout << bestHaps[i] << endl;
        for (size_t i = 0; i < 10; i++)
            cout << "0x" << hex << bestHapsFlags.getData()[i] << dec << endl;
        cout << endl;
    }
#endif
#endif

#ifdef STOPWATCH
    Stopwatch swcl("CallLocs/Monos");
#endif
    findMonos();
    findCallLocs();
#ifdef STOPWATCH
    swcl.stop();
#endif

#ifdef DEBUG_TARGET
#ifdef DEBUG_PRINT_LOCATIONS
    vector<size_t> monhets;
    for (auto call : callLocs) {
        if (monos[splitsites[call]])
            monhets.push_back(splitsites[call]);
    }
    cout << "Monos (#" << monhets.size() << "):" << endl;
    for (auto m : monhets)
        cout << " " << m;
    cout << "\nSplits:" << endl;
    for (auto split : splitsites)
        cout << " " << split;
    cout << "\nCalls:" << endl;
    for (auto call : callLocs)
        cout << " " << splitsites[call];
    cout << "\nHoms:" << endl;
    int cli = 0;
    for (auto s : splitsites) {
        if (s != splitsites[callLocs[cli]])
            cout << " " << s;
        else
            cli++;
    }
    cout << endl;
#endif
#endif
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
    cout << "target: " << ntarget << " splits: " << splitsites.size() << " (calls: " << callLocs.size() << " homs: " << (splitsites.size() - callLocs.size()) << ")" << endl;
#endif
}


void Target::phase() {

    // shortcut if target is haploid
    // NOTE: phasedTargetPat and phasedDosagePat are not required
    if (vcfdata.getHaploidsTgt()[ntarget]) {
    // copy target data to phasedTarget
        phasedDosageMat.clear();
        phasedDosageMat.resize(vcfdata.getNSNPs(), 0.0);
        for (size_t m = 0; m < MFull; m++) {
            Haplotype h = Haplotype::Ref; // default
            switch (targetFull[m]) {
            case Genotype::HomRef:
                break;
            case Genotype::HomAlt:
                h = Haplotype::Alt;
                phasedDosageMat[m] = 1.0;
                break;
            case Genotype::Miss:
                // "impute" a missing genotype according to RefPanelAF if "impMissing" is set, otherwise everything's fine since the output writer will mark this site as missing anyway
                if (impMissing) {
                    float af = vcfdata.getAlleleFreqCommon(m);
                    phasedDosageMat[m] = af;
                    h = af > 0.5 ? Haplotype::Alt : Haplotype::Ref;
                    totalconf += af > 0.5 ? af : 1.0-af;
                }
                totalmissings++;
                break;
            default: // heterozygous is not possible here
                cout << endl;
                StatusFile::addError("Found heterozygous genotype in haploid target.");
                exit(EXIT_FAILURE);
            }
            phasedTargetMat.push_back_withPreInit(toBool(h));
        }
        totalconf += MFull - totalmissings;
        ncallsites += MFull - totalmissings;
        return;
    } // end haploid target

    // shortcut if phasing should be skipped
    if (skipPhasing) {
    // copy target data to phasedTarget, use phase information from input file
        phasedDosageMat.clear();
        phasedDosageMat.resize(vcfdata.getNSNPs(), 0.0);
        phasedDosagePat.clear();
        phasedDosagePat.resize(vcfdata.getNSNPs(), 0.0);
        for (size_t m = 0; m < MFull; m++) {
            Haplotype hmat = Haplotype::Ref; // default
            Haplotype hpat = Haplotype::Ref; // default
            switch (targetFull[m]) {
            case Genotype::HomRef:
                break;
            case Genotype::HomAlt:
                hmat = hpat = Haplotype::Alt;
                phasedDosageMat[m] = 1.0;
                phasedDosagePat[m] = 1.0;
                break;
            case Genotype::Miss:
                // "impute" a missing genotype according to RefPanelAF if "impMissing" is set, otherwise everything's fine since the output writer will mark this site as missing anyway
                if (impMissing) {
                    float af = vcfdata.getAlleleFreqCommon(m);
                    phasedDosageMat[m] = af;
                    phasedDosagePat[m] = af;
                    hmat = hpat = af > 0.5 ? Haplotype::Alt : Haplotype::Ref;
                    totalconf += af > 0.5 ? af : 1.0-af;
                }
                totalmissings++;
                break;
            default: // heterozygous - phase according to input
                if (vcfdata.getTgtInPhases()[ntarget][m]) {
                    hmat = Haplotype::Alt;
                    phasedDosageMat[m] = 1.0;
                } else {
                    hpat = Haplotype::Alt;
                    phasedDosagePat[m] = 1.0;
                }
            }
            phasedTargetMat.push_back_withPreInit(toBool(hmat));
            phasedTargetPat.push_back_withPreInit(toBool(hpat));
        }
        totalconf += MFull - totalmissings;
        ncallsites += MFull - totalmissings;
        return;
    } // end skipPhasing

    prepPhase();

    size_t M = 2*splitsites.size()+1;
    size_t capacity = roundToMultiple(K, UNITWORDS * sizeof(BooleanVector::data_type) * 8) / 8;
    BooleanVector::data_type *refincdata = (BooleanVector::data_type*) MyMalloc::calloc(M * capacity, 1, string("refincdata_t")+to_string(ntarget)); // init complete area with false

    vector<BooleanVector> refincT(M, BooleanVector());
    auto curr_data = refincdata;
    for (auto &ref : refincT) {
        ref.setData(curr_data, capacity, K); // init with size K
        curr_data += capacity / sizeof(BooleanVector::data_type); // jump to next ref/incon (alternating)
    }

    vector<int>* count0 = new vector<int>(M);
    // make transposed condensed reference
#ifdef STOPWATCH
    Stopwatch swcref("Cref+Inc");
#endif
    condensedReferenceT(refincT, count0);
#ifdef STOPWATCH
    swcref.stop();
#endif
//
//#ifdef STOPWATCH
//    Stopwatch swinc("Incons+IBD");
//#endif
//    findInconsAndCheckIBD(refincT);
//#ifdef STOPWATCH
//    swinc.stop();
//#endif

    PBWT pbwt(&refincT, count0, ntarget);

    phaseExtPBWT(pbwt);

    MyMalloc::free(refincdata);

}


void Target::phaseExtPBWT(PBWT &pbwt) {

#ifdef STOPWATCH
    Stopwatch phaseextsw("PhaseExtPBWT");
#endif

    const auto &csplits(splitsites);

    // create the constraints vector for path extension and fill unconstrained for each split site
    // a constraint targets the site itself in the homozygous case,
    // for the heterozygous case, it addresses the segment _before_ the site ('left' in forward order)
    vector<Constraint> constraints;
    constraints.resize(csplits.size(), Constraint(Constraint::Type::NoConstr, 0));

    // resize the dosage storage
    phasedDosageMat.clear();
    phasedDosagePat.clear();
    phasedDosageMat.resize(MFull, 0.0);
    phasedDosagePat.resize(MFull, 0.0);

    Haplotype currMatPhase = vcfdata.getStartPhase(ntarget); // set the phase of the first het as indicated from the last chunk
    Haplotype lastMatPhase = vcfdata.getStartPhase(ntarget); // also need to store the previous phase which we set to the same phase w.l.o.g.

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
    vector<fp_type> fastProbs;
    fastProbs.reserve(callLocs.size()); // probs for segment before each callloc
    fastProbs.push_back(0); // first decision is fixed to 0/1, but prob is fifty-fifty (shifted by -0.5)
    vector<Haplotype> fastHaps;
    fastHaps.reserve(callLocs.size());
    fastHaps.push_back(currMatPhase);
#endif


    if (doPrePhasing) { // fast pre-phasing
#ifdef STOPWATCH
        Stopwatch sw("Phase Fast");
#endif

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
        cout << "\nFAST PRE-PHASING\n" << endl;
#endif

        size_t history = static_cast<size_t>((HistorySizeFast * hist) + 0.5);
        Beam beam(history, BeamWidthFast, DeltaFast, //DeltaFast+2,
                pErr, pLimit, expectIBDcM,
                cMsplitDiffs, callLocs, target, pbwt, constraints);

        // the first call site gets a 'heterozygous keep' constraint such that only hets will be created at the very first site
        if (csplits.size() > 0)
            constraints[0] = Constraint(Constraint::Type::HetKeep, 0);
        MinMaxHeap<ConstraintContainer> callheap;
        MinMaxHeap<ConstraintContainer> homheap;
        // starting at second call loc (first will be handled automatically by advancing the beam front and is het 0/1 anyway)
        size_t next_call_idx = 1;
        for (size_t si = 1; si < csplits.size(); si++) { // over all split sites, starting with second
            if (si == callLocs[next_call_idx]) { // split site is a call site, so call probability and insert into constraint heap

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                cout << "\ni: " << next_call_idx << " s: " << csplits[si] << " phase: " << currMatPhase << " target: " << target[si] << " constr: " << constraints[si] << endl;
#endif
                fp_type p = beam.callKeepPhaseProb(next_call_idx) - 0.5; // shift value in area [-0.5:0.5]
                Haplotype ref = currMatPhase; // reference haplotype is the current phase
                Constraint c(Constraint::Type::NoConstr, callLocs[next_call_idx] - callLocs[next_call_idx-1] - 1); // the type is decided later, but the reference is set to one call loc back per default (which is 0 if there are no hom splits in between)
                if (next_call_idx > 1) {
                    fp_type p2 = beam.callKeepPhaseProb(next_call_idx, true) - 0.5; // shift value in area [-0.5:0.5]
                    if (abs(p2) > abs(p)) { // which one is more distant to 50%?
                        p = p2; // take the probability which is more definite
                        ref = lastMatPhase; // reference whether to switch or keep is two sites back!
                        c.refidx = callLocs[next_call_idx] - callLocs[next_call_idx-2] - 1; // correct the reference index for the constraint
                    }
                }

                lastMatPhase = currMatPhase; // save old phase
                if (p < 0) { // switch phase according to reference
                    currMatPhase = !ref;
                    c.type = Constraint::Type::HetSwitch;
                } else { // keep phase according to reference
                    currMatPhase = ref;
                    c.type = Constraint::Type::HetKeep;
                }
                constraints[si] = c; // replace constraint
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                fastProbs.push_back(p);
                fastHaps.push_back(currMatPhase);
                cout << " p: " << (p+0.5) << " new phase: " << currMatPhase << " (" << (p<0 ? "switch" : "keep") << ") new constr: " << constraints[si] << endl;
#endif

                // insert constraint into heap
                p = abs(p); // only positives for the heap (is now a shifted confidence)
                callheap.emplace(si, p);

                next_call_idx++;


            } else { // homozygous split site -> constraint for hom heap

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                cout << "\nhom split: " << csplits[si] << " target: " << target[si] << " constr: " << constraints[si] << endl;
#endif
                fp_type p = beam.callDosage(si); // diploid dosage in [0.0:2.0]
                if (target[si] == Genotype::HomRef) {
                    constraints[si] = Constraint(Constraint::Type::HomRef, 0); // replace constraint
                    p = 2.0 - p; // the dosage confidence is inverse to the dosage: a dosage of 0.0 indicates a high confidence (2.0) for a homozygous reference genotype, and the other way round
                } else // has to be Genotype::HomAlt -> don't need to adjust dosage, already corresponds to dosage confidence
                    constraints[si] = Constraint(Constraint::Type::HomAlt, 0); // replace constraint
                homheap.emplace(si, p);
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                cout << " tgt.dosage: " << p << " new constr: " << constraints[si] << endl;
#endif
            }
        } // END over all split sites

        // The heaps now contain all constraints sorted by their confidence.
        // According to Eagle at least 50% will be fixed. If more than 50% have a confidence higher than 99%, those will be fixed as well.
        size_t sizedivdmin = callheap.getSize()*minFixSize;
        size_t sizedivdmax = callheap.getSize()*(1.0-maxFixSize);
        size_t si;
        for (si = 0; si < sizedivdmin; si++) { // iterate over the minimum half of the heap to reset those constraints
            ConstraintContainer cc = callheap.popMin();
            if (cc.p >= maxFixThresh-0.5 || (cc.p >= minFixThresh-0.5 && si >= sizedivdmax)) // -0.5 to make the threshold comparable with the shifted confidences (in [0.0:0.5]) as stored in the constraint containers
                break; // all remaining constraints are more confident and stay fixed
            // else: reset constraint because it is in lower half and below 99% or greater 99% but not 100% and in lowest 10%
            constraints[cc.constraint_idx].type = Constraint::Type::NoConstr;
        }
        if (si == sizedivdmin) { // all constraints until here were removed
            // go on removing fifty-fifty constraints
            for (; si < callheap.getSize(); si++) {
                ConstraintContainer cc = callheap.popMin();
                if (cc.p > 0.0) // -0.5 to make the threshold comparable with the shifted confidences (in [0.0:0.5]) as stored in the constraint containers
                    break; // all remaining constraints are more confident and stay fixed
                // else: reset constraint because it's fifty-fifty
                constraints[cc.constraint_idx].type = Constraint::Type::NoConstr;
            }
        }

        // the same for the homozygous sites
        sizedivdmin = homheap.getSize()*minFixSize;
        sizedivdmax = homheap.getSize()*(1.0-maxFixSize);
        for (si = 0; si < sizedivdmin; si++) { // iterate over the minimum half of the heap to reset those constraints
            ConstraintContainer cc = homheap.popMin();
            if (cc.p >= maxFixThresh*2.0 || (cc.p >= minFixThresh*2.0 && si >= sizedivdmax)) // *2.0 to make the threshold comparable with diploid dosages (in [0.0:2.0]) as stored in the constraint containers
                break; // all remaining constraints are more confident and stay fixed
            // else: reset constraint because it is in lower half and below 99% or greater 99% but not 100% and in lowest 10%
            constraints[cc.constraint_idx].type = Constraint::Type::NoConstr;
        }
        // no need to further remove fifty-fifty's

#ifdef STOPWATCH
        sw.stop();
#endif
    } // END fast pre-phasing

    currMatPhase = vcfdata.getStartPhase(ntarget); // start over new
    // store the per-segment phase probabilities for later comparison in the reverse post-phasing step
    vector<fp_type> phaseProbs;
    phaseProbs.reserve(callLocs.size()); // probs for segment before each callloc
    phaseProbs.push_back(0); // first decision is fixed to 0/1, but prob is fifty-fifty (shifted by -0.5)
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
    vector<Haplotype> fwdHaps;
    fwdHaps.reserve(callLocs.size());
    fwdHaps.push_back(currMatPhase);
    int misscount = 0;
    int monhetcount = 0;
#endif
    {
#ifdef STOPWATCH
        Stopwatch sw("Phase Fine + Call");
#endif

        // constraints are set now, so we can run the fine search
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
        cout << "\nFINE PHASING AND PHASE CALLING\n" << endl;
#endif

        // do the phase calls according to the probabilities from beam.callKeepPhaseProb()

        size_t history = static_cast<size_t>((HistorySize * hist) + 0.5);
        Beam beam(history, BeamWidth, Delta, //Delta+1,
                pErr, pLimit, expectIBDcM,
                cMsplitDiffs, callLocs, target, pbwt, constraints);

        // pre-processing: set all haplotypes according to homozygous genotypes until first het
        size_t m = 0; // absolute site cursor
        while (csplits.size() > 0 ? m < csplits[callLocs[0]] : m < MFull) { // until first heterozygous call location (callLocs[0] should always be zero, so one could also write splits[0])
            // default is homozygous reference
            Haplotype mat = Haplotype::Ref;
            Haplotype pat = Haplotype::Ref;
            switch (targetFull[m]) {
            case Genotype::HomAlt:
                mat = Haplotype::Alt;
                pat = Haplotype::Alt;
                phasedDosageMat[m] = 1.0;
                phasedDosagePat[m] = 1.0;
                break;
            case Genotype::Miss:
                // "impute" a missing genotype according to RefPanelAF if "impMissing" is set, otherwise everything's fine since the output writer will mark this site as missing anyway
                if (impMissing) {
                    float af = vcfdata.getAlleleFreqCommon(m);
                    phasedDosageMat[m] = af;
                    phasedDosagePat[m] = af;
                    mat = af > 0.5 ? Haplotype::Alt : Haplotype::Ref;
                    pat = mat;
                }
                if (lastIter) { // only count in last iteration
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                    misscount++;
#endif
                    totalmissings++;
                }
                break;
            case Genotype::Het: // heterozygous -> must be mono, otherwise it would have been a call location
                mat = Haplotype::Ref;
                pat = Haplotype::Alt;
                phasedDosageMat[m] = 0.5;
                phasedDosagePat[m] = 0.5;
                break;
            default: // nothing to do for hom ref
                break;
            }
            phasedTargetMat.push_back_withPreInit(toBool(mat));
            phasedTargetPat.push_back_withPreInit(toBool(pat));
            m++;
        }

        // init first het
        if (m < MFull) { // there are hets in the target (usual case, but the other case is also possible)
            phasedTargetMat.push_back_withPreInit(toBool(currMatPhase));
            phasedTargetPat.push_back_withPreInit(!toBool(currMatPhase));
            // dosage will be set by reverse run
            m++; // now points to site after first het
        }

        size_t next_call_idx = 1;
        // iterate over all remaining sites
        for (; m < MFull; m++) {
            // default is homozygous reference
            Haplotype mat = Haplotype::Ref;
            Haplotype pat = Haplotype::Ref;
            switch (targetFull[m]) {
            case Genotype::Het: // heterozygous, choose the haplotype according to the keep-phase probability
                if (!monos[m]) { // call phase only for non-mono sites
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                    cout << "\ni: " << next_call_idx << " s: " << m << " phase: " << currMatPhase << " target: " << target[callLocs[next_call_idx]] << " constr: " << constraints[callLocs[next_call_idx]] << endl;
#endif
                    fp_type p = beam.callKeepPhaseProb(next_call_idx) - 0.5; // shift value in area [-0.5:0.5]
                    Haplotype ref = currMatPhase; // reference haplotype is the current phase

                    // store phase probability for segment
                    phaseProbs.push_back(p);
                    if (p < 0) // switch phase according to reference
                        currMatPhase = !ref;
                    else // keep phase according to reference
                        currMatPhase = ref;
                    next_call_idx++;

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                    fwdHaps.push_back(currMatPhase);
#endif
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                    cout << " p: " << (p+0.5) << " new phase: " << currMatPhase << " (" << (p<0 ? "switch" : "keep") << ")" << endl;
#endif
                } else { // mono site: keep phase, but set dosage to 0.5
                    phasedDosageMat[m] = 0.5;
                    phasedDosagePat[m] = 0.5;
                    if (lastIter) { // only count in last iteration
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                        monhetcount++;
#endif
                        totalmonhets++;
                    }

                }
                mat = currMatPhase;
                pat = !currMatPhase;
                break;
            case Genotype::HomAlt:
                mat = Haplotype::Alt;
                pat = Haplotype::Alt;
                phasedDosageMat[m] = 1.0;
                phasedDosagePat[m] = 1.0;
                break;
            case Genotype::Miss:
                // "impute" a missing genotype according to RefPanelAF if "impMissing" is set, otherwise everything's fine since the output writer will mark this site as missing anyway
                if (impMissing) {
                    float af = vcfdata.getAlleleFreqCommon(m);
                    phasedDosageMat[m] = af;
                    phasedDosagePat[m] = af;
                    mat = af > 0.5 ? Haplotype::Alt : Haplotype::Ref;
                    pat = mat;
                }
                if (lastIter) { // only count in last iteration
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                    misscount++;
#endif
                    totalmissings++;
                }
                break;
            default: // nothing to do for hom ref
                break;
            }
            phasedTargetMat.push_back_withPreInit(toBool(mat));
            phasedTargetPat.push_back_withPreInit(toBool(pat));
        }
#ifdef STOPWATCH
        sw.stop();
#endif
    }

    // reverse fine phasing (only in the last (or only) iteration, if desired)
    // the current phase from the end of the last forward run is kept and used as the init phase for the reverse run

    // store the per-call-site confidences (in [0.5,1.0])
    vector<fp_type> callConfs;
    if (lastIter && csplits.size() > 0) {
        callConfs.resize(callLocs.size(), 0.0);
        if (doRevPhasing) {
            // the last phase has been decided by the call for the last segment in the fwd run,
            // all others will be re-calculated and compared to the previous one before being stored during reverese phasing
            callConfs[callLocs.size()-1] = abs(phaseProbs.back()) + 0.5;
        } else {
            // no reverse phasing -> copy all previous phase probabilities
            auto cc = callConfs.begin();
            for (const auto &p : phaseProbs) {
                *cc = abs(p)+0.5;
                cc++;
            }
        }
    }
    // reverse phasing
    if (lastIter && doRevPhasing && csplits.size() > 0) {
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
        vector<fp_type> revProbs(callLocs.size()); // probs for segment before each callloc (here in reverse order)
        revProbs[callLocs.size()-1] = phaseProbs.back(); // the first reverse prob is simply taken over by the forward run (doesn't matter anyway)
        vector<Haplotype> revHaps(callLocs.size());
        revHaps[callLocs.size()-1] = currMatPhase; // take over the current phase from the forward run
#endif

        // the phase of the hets is improved by a reverse post-phasing step

#ifdef STOPWATCH
        Stopwatch swrpbwt("revPBWT");
#endif
        // set PBWT to reverse
        pbwt.switchToReverse();
#ifdef STOPWATCH
        swrpbwt.stop();
#endif

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
        cout << "\nREVERSE FINE POST-PHASING\n" << endl;
#endif
#ifdef STOPWATCH
        Stopwatch sw("Phase Rev");
#endif
        size_t history = static_cast<size_t>((HistorySize * hist) + 0.5);

        // reverse constraints:
        // a constraint targets the site itself in the homozygous case,
        // for the heterozygous case, it addresses the segment _before_ the site ('right' in reverse order)
        vector<Constraint> revconstraints;
        revconstraints.resize(csplits.size(), Constraint(Constraint::Type::NoConstr, 0));
        for (size_t cidx = 1; cidx < csplits.size(); cidx++) { // starting from second (first is the init constraint not required by the reverse run)
            Constraint &c = constraints[cidx];
            switch (c.type) {
            case Constraint::Type::HetKeep:
            case Constraint::Type::HetSwitch:
                {
                    Constraint &rc = revconstraints[cidx-c.refidx-1]; // where the original constraint refers to
                    if (rc.type != Constraint::Type::NoConstr) { // there was a constraint already set (must be het!)
                        // follow the constraint and apply another constraint refering to this site
                        size_t ridx = cidx - c.refidx + rc.refidx;
                        revconstraints[ridx].refidx = c.refidx - rc.refidx - 1;
                        if (c.type != rc.type) // need to set as switch
                            revconstraints[ridx].type = Constraint::Type::HetSwitch;
                        else // need to set as keep
                            revconstraints[ridx].type = Constraint::Type::HetKeep;
                    } else { // copy the constraint
                        //rc = c; // I don't trust this if it is really copying...
                        rc.type = c.type;
                        rc.refidx = c.refidx;
                    }
                }
                break;

            case Constraint::Type::HomRef:
            case Constraint::Type::HomAlt:
                // simply copy constraint type, no adaption of reference index required
                revconstraints[cidx].type = c.type;
                break;
            default: // nothing to do for unconstrained
                break;
            }
        }
        revconstraints[csplits.size()-1] = Constraint(Constraint::Type::HetKeep, 0); // last constraint will not be set by above procedure and is the init constraint for the reverse run

        // store dosage for last het (reverse run recalculates only from the second last)
        float dref = phaseProbs.back() < 0.0 ? 0.5 + phaseProbs.back() : 0.5 - phaseProbs.back(); // dosage for reference allele (between 0.0 and 0.5)
        float dalt = phaseProbs.back() < 0.0 ? 0.5 - phaseProbs.back() : 0.5 + phaseProbs.back(); // dosage for alternative allele (between 0.5 and 1.0)
        phasedDosageMat[csplits[callLocs.back()]] = currMatPhase == Haplotype::Ref ? dref : dalt;
        phasedDosagePat[csplits[callLocs.back()]] = currMatPhase == Haplotype::Ref ? dalt : dref;

        // Here the mat phase of the beam is w.l.o.g. initialized with Haplotype::Ref, although the last "currMatPhase" from the forward run could be Haplotype::Alt.
        // This is ok since both hap paths will behave the same way, no matter if they are maternal or paternal.
        Beam beam(history, BeamWidth, Delta, //Delta+1,
                pErr, pLimit, expectIBDcM,
                cMsplitDiffs, callLocs, target, pbwt, revconstraints,
                true); // reverse Beam!!

        // starting at second reverse call loc (first will be handled automatically by advancing the beam front and is het 0/1 anyway)
        for (int rci = callLocs.size()-2; rci >= 0; rci--) { // over all call sites, starting with second (in reverse order)
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
            cout << "\ni: " << rci << " s: " << csplits[callLocs[rci]] << " phase: " << currMatPhase << " target: " << target[callLocs[rci]] << " constr: " << revconstraints[callLocs[rci]] << endl;
#endif
            fp_type p = beam.callKeepPhaseProb(rci) - 0.5; // shift value in area [-0.5:0.5]
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
            revProbs[rci] = p;
#endif
            Haplotype ref = currMatPhase; // reference haplotype is the current phase
            // test for the better estimation from fwd and rev run
            // ATTENTION on indices! probs are stored for the segment _before_ the call loc (i.e. left for forward, right for backward)
            if (abs(phaseProbs[rci+1]) > abs(p))
                p = phaseProbs[rci+1];
            if (p < 0) { // switch phase according to reference
                currMatPhase = !ref;
            } else { // keep phase according to reference
                currMatPhase = ref;
            }
            // store phase
            phasedTargetMat.set(csplits[callLocs[rci]], toBool(currMatPhase));
            phasedTargetPat.set(csplits[callLocs[rci]], !toBool(currMatPhase));
            // store dosage
            dref = p < 0.0 ? 0.5 + p : 0.5 - p; // dosage for reference allele (between 0.0 and 0.5)
            dalt = p < 0.0 ? 0.5 - p : 0.5 + p; // dosage for alternative allele (between 0.5 and 1.0)
            phasedDosageMat[csplits[callLocs[rci]]] = currMatPhase == Haplotype::Ref ? dref : dalt;
            phasedDosagePat[csplits[callLocs[rci]]] = currMatPhase == Haplotype::Ref ? dalt : dref;
            // store confidence (shifted back to [0.5:1])
            callConfs[rci] = abs(p) + 0.5;
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
            revHaps[rci] = currMatPhase;
#endif
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
            cout << " p: " << (p+0.5) << " new phase: " << currMatPhase << " (" << (p<0 ? "switch" : "keep") << ")" << endl;
#endif

        } // END over all call sites
#ifdef STOPWATCH
        sw.stop();
#endif

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
        cout << "\nPer-site phase probabilities for target " << ntarget << ":" << endl;
        cout << "site\tfst.prob\tfwd.prob.\trev.prob.\tfst.bit\tfwd.bit\trev.bit\tf.con.\tr.con\tdiffAB\tdiffBC" << endl;

        size_t ci;
        for (ci = 0; ci < callLocs.size(); ci++) {
            cout << csplits[callLocs[ci]] << "\t";
            cout << setprecision(6) << setw(8) << fastProbs[ci]+0.5 << "\t"; // fast phase prob from this site to next
            cout << setprecision(6) << setw(8) << phaseProbs[ci]+0.5 << "\t"; // phase prob from this site to next
            cout << setprecision(6) << setw(8) << revProbs[ci]+0.5 << "\t";   // reverse phase prob from next site to this
            cout << (fastHaps[ci] == Haplotype::Alt ? "1" : "0") << "\t";
            cout << (fwdHaps[ci] == Haplotype::Alt ? "1" : "0") << "\t";
            cout << (revHaps[ci] == Haplotype::Alt ? "1" : "0") << "\t";
            cout << constraints[callLocs[ci]] << "\t" << revconstraints[callLocs[ci]] << "\t";
            cout << (fastHaps[ci] != fwdHaps[ci] ? "X" : " ") << "\t";
            cout << (fwdHaps[ci] != revHaps[ci] ? "X" : " ") << endl;
        }
#endif
    } // END reverse fine phasing

    // calculate total confidence for sample
    if (lastIter && csplits.size() > 0) {
        for (size_t cci = 0; cci < callConfs.size(); cci++) {
            totalconf += callConfs[cci];
        }
        ncallsites += callConfs.size();

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
        cout << "Missings: " << misscount << " Heterozygous singletons: " << monhetcount << endl;
        cout << "Average phase confidence for target " << ntarget << ": " << totalconf/ncallsites << endl;
#endif
    }

#ifdef STOPWATCH
    phaseextsw.stop();
#endif

}

void Target::findBestHaps() {
    bestHaps.clear();
    bestHaps.resize(K);
//    // dummy element (required e.g. for IBD check)
//    bestHaps[K] = 0ull;

    if (useExtArch) {
        size_t bhcapacity = roundToMultiple(NrefhapsAllDiploid, UNITWORDS * sizeof(BooleanVector::data_type) * 8) / 8;
        BooleanVector::data_type *bhdata = (BooleanVector::data_type*) MyMalloc::calloc(bhcapacity, 1, string("bhdata_t")+to_string(ntarget));
        bestHapsFlags.setData(bhdata, bhcapacity, NrefhapsAllDiploid);
    }

    // if we are in an iteration > 1, this is the index of this target in the reference, which has to be ignored
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
    size_t ignoreidx = ntarget-DEBUG_TARGET_START + vcfdata.getNReferenceHaps()/2;
#else
    size_t ignoreidx = ntarget + vcfdata.getNReferenceHaps()/2;
#endif

    if (K >= NrefhapsCorrected) { // simply take all references
        if (NrefhapsAllDiploid != NrefhapsCorrected) { // need to jump over haps and/or ignore target
            for (size_t n = 0, nr = 0; nr < NrefhapsAllDiploid/2; nr++) {
                if (nr != ignoreidx) {
                    bestHaps[n++] = 2*nr;
                    if (!vcfdata.getHaploidsRef()[nr])
                        bestHaps[n++] = 2*nr+1;
                }
            }
            if (useExtArch) { // create best haps flags
                for (auto bh : bestHaps)
                    bestHapsFlags.setWithPreInit(bh, true);
            }
        } else { // neither haps nor a target to ignore
            for (size_t nr = 0; nr < NrefhapsAllDiploid; nr++) {
                bestHaps[nr] = nr;
            }
            if (useExtArch) { // create best haps flags
                // simply set all NrefhapsAllDiploid flags to 1
                memset(bestHapsFlags.getData(), 0xff, (NrefhapsAllDiploid/(sizeof(BooleanVector::data_type)*8))*sizeof(BooleanVector::data_type));
                // remainder
                size_t remainder = NrefhapsAllDiploid%(sizeof(BooleanVector::data_type)*8);
                if (remainder) {
                    for (size_t nr = (NrefhapsAllDiploid/(sizeof(BooleanVector::data_type)*8))*(sizeof(BooleanVector::data_type)*8); nr < NrefhapsAllDiploid; nr++) {
                        bestHapsFlags.setWithPreInit(nr, true);
                    }
                }
            }
        }
    } else { // K < num haps
        vector<pair<size_t,size_t>> errors;
        errors.reserve(NrefhapsCorrected);
        for (size_t nr = 0; nr < NrefhapsAllDiploid/2; nr++) {
            if (nr != ignoreidx)
            {
                size_t nerr = targetFull.compare_withPreInit(referenceFull[2*nr]);
                errors.push_back(make_pair(nerr, 2*nr));
                if (!vcfdata.getHaploidsRef()[nr]) {
                    nerr = targetFull.compare_withPreInit(referenceFull[2*nr+1]);
                    errors.push_back(make_pair(nerr, 2*nr+1));
                }
            }
        }
        sort(errors.begin(), errors.end());
        for (size_t k = 0; k < K; k++)
            bestHaps[k] = errors[k].second;

        // sorting by index is only required for IBD check, anyway it is more memory friendly and costs almost nothing
//        sort(bestHaps.begin(), bestHaps.end()-1); // sort without the dummy! (required for IBD check)
        sort(bestHaps.begin(), bestHaps.end());

        // set flags
        if (useExtArch) {
            for (auto bh : bestHaps)
                bestHapsFlags.setWithPreInit(bh, true);
        }

//        // DEBUG
//        vector<pair<size_t,size_t>> bestHapsPairs;
//        vector<pair<size_t,size_t>> notInList;
//        if (ntarget == 0) {
//            for (size_t k = 0; k < K; k++)
//                bestHapsPairs.push_back(make_pair(errors[k].second, errors[k].first));
//            for (size_t k = K; k < 2*NrefhapsCorrected; k++)
//                notInList.push_back(make_pair(errors[k].second, errors[k].first));
//            sort(bestHapsPairs.begin(), bestHapsPairs.end());
//            sort(notInList.begin(), notInList.end());
//
//            cout << "best haps: " << endl;
//            for (size_t k = 0; k < K; k++)
//                cout << bestHapsPairs[k].first << ": " << bestHapsPairs[k].second << endl;
//            cout << "not in list: " << endl;
//            for (size_t k = K; k < 2*NrefhapsCorrected; k++)
//                cout << notInList[k-K].first << ": " << notInList[k-K].second << endl;
//        }
//        // __DEBUG
    }
}


void Target::findMonos() {
    size_t capacity = roundToMultiple(MFull, UNITWORDS * sizeof(BooleanVector::data_type) * 8) / 8;
    BooleanVector::data_type *monodata = (BooleanVector::data_type*) MyMalloc::malloc(capacity, string("monodata_t")+to_string(ntarget));
    monos.setDataAndInit(monodata, capacity, MFull, true); // initialize all sites as mono
    for (size_t curr = 1; curr < K; curr++) { // over all referencing haplotypes
        const auto &lastdata = referenceFull[bestHaps[curr-1]].getData();
        const auto &currdata = referenceFull[bestHaps[curr]].getData();
        size_t limit = capacity / sizeof(BooleanVector::data_type);
        for (size_t i = 0; i < limit; i++) {
            uint64_t compare = ~(lastdata[i] ^ currdata[i]); // equal positions are now '1'
            monodata[i] &= compare; // mono sites until here are still marked with '1'
        }
    }
//    // DEBUG
//    cerr << "Monos: " << MFull << endl;
//    for (size_t i = 0; i < MFull; i++) {
//        cerr << (monos[i] ? "1" : "0");
//    }
//    cerr << endl;
}


void Target::findCallLocs() {

    size_t capacity = roundToMultiple(MFull, UNITWORDS * sizeof(BooleanVector::data_type) * 8) / 8;
    BooleanVector::data_type *splitdata = (BooleanVector::data_type*) MyMalloc::malloc(capacity, string("splitdata")+to_string(ntarget));
    splits.setDataAndInit(splitdata, capacity, MFull, false); // default: all sites are not split sites

    fp_type lastSplitcM = cMpos[0];
    // first segment to phase is always from first het, so find first het and go on from there
    size_t m;
    for (m = 0; m < MFull; m++) {
        //if (!monos[m] && targetFull[m] == Genotype::Het)
        if (targetFull[m] == Genotype::Het) // since the first het site will be set to 0,1 w.l.o.g., it may also be a mono site.
            break;
    }
    for (; m < MFull; m++) {
        if (!monos[m] && targetFull[m] != Genotype::Miss // not mono and not missing
                && (targetFull[m] == Genotype::Het || cMpos[m] - lastSplitcM > cMmaxSplit)) { // and heterozygous or last split is too far away
            if(targetFull[m] == Genotype::Het) // heterozygous -> call location!
                callLocs.push_back(splitsites.size()); // i.e. the index in splits where this site will be found after the next command
            splitsites.push_back(m);
            target.push_back(targetFull[m]);
            splits.setWithPreInit(m, true); // mark as split site!
            lastSplitcM = cMpos[m];
        }
    }

    // simply return if we don't have a split site (i.e. without or only mono hets)
    if (splitsites.size() == 0)
        return;

    // also return, if there are no call sites (but cleanup split sites that may result from splitting long segments)
    if (callLocs.size() == 0) {
        splitsites.clear();
        splits.setDataAndInit(splitdata, capacity, MFull, false); // easiest way: simply reset the whole structure
        return;
    }

    // truncate splits such that the last split site is a (heterozygous) call site,
    // all trailing homozygous sites are erased
    splitsites.resize(callLocs.back()+1);
    splits.clearAfter(splitsites.back());

    // fill split site differences
    const auto &csplits(splitsites);
    cMsplitDiffs.reserve(csplits.size()+1);
    cMsplitDiffs.push_back(cMpos[csplits[0]] - cMpos[0]); // first diff is the length of the segment to the first split site
    for (size_t i = 1; i < csplits.size(); i++) // starting with '1': diff[1] = split[1]-split[0];
        cMsplitDiffs.push_back(cMpos[csplits[i]] - cMpos[csplits[i-1]]);
    cMsplitDiffs.push_back(cMpos.back() - cMpos[csplits.back()]); // last diff is the length of the segment after the last split site

}

// this is the second most computationally expensive part of the phasing!
void Target::condensedReferenceT(vector<BooleanVector> &refincT, vector<int> *count0) {

    const auto &csplits(splitsites);

    size_t inconFullTcap = roundToMultiple(NrefhapsAllDiploid, UNITWORDS * sizeof(BooleanVector::data_type) * 8) / 8;
    BooleanVector::data_type* inconFulldataT = (BooleanVector::data_type*) MyMalloc::calloc(inconFullTcap, 1, string("inconFulldataT_t")+to_string(ntarget));
    BooleanVector inconFullT(inconFulldataT, inconFullTcap, NrefhapsAllDiploid);

    auto currrefincTit = refincT.begin();
    auto csplitit = csplits.cbegin();
    size_t nextsplit = csplitit == csplits.cend() ? ~0ull : *csplitit++;

    size_t nrefhapsADwordsnopad = (NrefhapsAllDiploid+8*sizeof(BooleanVector::data_type)-1)/(8*sizeof(BooleanVector::data_type));
    BooleanVector::data_type nrefhapsADpadmask = (NrefhapsAllDiploid % (8*sizeof(BooleanVector::data_type))) ? ((1ull << (NrefhapsAllDiploid % (8*sizeof(BooleanVector::data_type)))) - 1) : ~(0ull);

    auto c0it = count0->begin();
    for (size_t m = 0; m < MFull; m++) { // for all variants

        if (m == nextsplit) { // reached split site
            // copy current incon data
            int count1 = 0;
            if (K == NrefhapsAllDiploid) { // copy all data directly
                // need to mask potentially false set padding bits
                inconFulldataT[nrefhapsADwordsnopad-1] &= nrefhapsADpadmask;
                memcpy(currrefincTit->getData(), inconFulldataT, nrefhapsADwordsnopad*sizeof(BooleanVector::data_type));
                for (size_t w = 0; w < nrefhapsADwordsnopad; w++)
                    count1 += __builtin_popcountll(inconFulldataT[w]);
            } else { // K is a subset
                for (size_t k = 0; k < K; k++) { // for all references
                    currrefincTit->setWithPreInit(k, inconFullT[bestHaps[k]]);
                    if (inconFullT[bestHaps[k]]) // inconsistent
                        count1++;
                }
            }
            *c0it = K - count1;
            c0it++;
            currrefincTit++; // points to ref (odd) part now

            // clean incon data for next segment
            memset(inconFulldataT, 0, inconFullTcap);

            // copy current reference data
            count1 = 0;
            if (K == NrefhapsAllDiploid) { // copy all data directly
                memcpy(currrefincTit->getData(), referenceFullT[m].getData(), nrefhapsADwordsnopad*sizeof(BooleanVector::data_type));
                for (size_t w = 0; w < nrefhapsADwordsnopad; w++)
                    count1 += __builtin_popcountll(referenceFullT[m].getData()[w]);
            } else { // K is a subset
                for (size_t k = 0; k < K; k++) { // for all references
                    currrefincTit->setWithPreInit(k, referenceFullT[m][bestHaps[k]]);
                    if (referenceFullT[m][bestHaps[k]]) // reference is 1 here
                        count1++;
                }
            }
            *c0it = K - count1;
            c0it++;
            currrefincTit++; // points to incon (even) part now

            // find next split site
            nextsplit = csplitit == csplits.cend() ? ~0ull : *csplitit++;

        } else { // site between splits
            // generate new incon data and combine with current
            mergeInconsistencies(targetFull[m], referenceFullT[m].getData(), inconFulldataT, nrefhapsADwordsnopad);
        }
    } // END for all variants

    // copy incon data from last segment
    int count1 = 0;
    if (K == NrefhapsAllDiploid) { // copy all data directly
        // need to mask potentially false set padding bits
        inconFulldataT[nrefhapsADwordsnopad-1] &= nrefhapsADpadmask;
        memcpy(currrefincTit->getData(), inconFulldataT, nrefhapsADwordsnopad*sizeof(BooleanVector::data_type));
        for (size_t w = 0; w < nrefhapsADwordsnopad; w++)
            count1 += __builtin_popcountll(inconFulldataT[w]);
    } else { // K is a subset
        for (size_t k = 0; k < K; k++) { // for all references
            currrefincTit->setWithPreInit(k, inconFullT[bestHaps[k]]);
            if (inconFullT[bestHaps[k]]) // inconsistent
                count1++;
        }
    }
    *c0it = K - count1;

    MyMalloc::free(inconFulldataT);
}

void Target::mergeInconsistencies(Genotype g, const BooleanVector::data_type* src, BooleanVector::data_type* dest, size_t size) {
    // determine inconsistency according to the current genotype
    switch (g) {
    case Genotype::HomRef:
        // all 1-haps are inconsistent to the 0-genotype
        for (size_t i = 0; i < size; i++) {
            dest[i] |= src[i];
        }
        break;
    case Genotype::HomAlt:
        // all 0-haps are inconsistent to the 1-genotype
        for (size_t i = 0; i < size; i++) {
            dest[i] |= ~(src[i]);
        }
        break;
    default:
        // missing or heterozygous genotypes are consistent to everything, so do nothing
        break;
    }
}

