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

#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <array>
#include <unordered_map>
#include <memory>

#include "MinMaxHeap.h"
#include "PBWT.h"
#include "BeamPath.h"
#include "StatusFile.h"

#include "Beam.h"

/* static */ const int Beam::minProbExp = -8;

bool operator<(const BeamPathContainer &a, const BeamPathContainer &b) { return a.prob < b.prob; }

Beam::Beam(size_t historySize, size_t beamWidth, size_t deltaCall,
        fp_type pErr_, fp_type pLimit_, fp_type expectIBDcM_,
        const vector<fp_type> &cMsplitDiffs_, const vector<size_t> &callLocs_,
        const GenotypeVector &target_,
        PBWT &pbwt_,
        const vector<Constraint> &constraints_,
        bool reverseBeam_)
    : H(historySize), P(beamWidth), DeltaCall(deltaCall),
      pErr(pErr_), pLimit(pLimit_), expectIBDcM(expectIBDcM_),
      cMsplitDiffs(cMsplitDiffs_), callLocs(callLocs_),
      target(target_),
      constraints(constraints_),
      pbwt(pbwt_),
      reverseBeam(reverseBeam_)
{
    prec.clear();
    prec.resize(historySize);

    cursor = -1; // before the first site
    curridx = -1;
    refidx = -1;
    if (reverseBeam)
        cursor = constraints.size(); // after the last site

    initBeamFront();
}

Beam::~Beam() {
    for (auto hp : happaths)
        delete hp;
}

void Beam::initBeamFront() {
    // calculate everything for first split site, i.e. initialize beam and set start probabilities

    // The first path is an empty path which is only filled pro-forma with a heterozygous site

    // create initial beam path with mat = Ref, pat = Alt
    HaplotypePath *mat = new HaplotypePath(H, pbwt, Haplotype::Ref);
    HaplotypePath *pat = new HaplotypePath(H, pbwt, Haplotype::Alt);
    happaths.clear();
    happaths.reserve(2*P); // at most there will be 2*beamWidth haplotype paths
    // correct haplotype order for happaths
    happaths.push_back(mat);
    happaths.push_back(pat);

    BeamPath bp(pErr, happaths, 0, 1, 1.0, 1.0, 0);
//    BeamPath bp(pErr, happaths, 1, 0, 1.0, 1.0, 0); // DEBUG: switch mat + pat for start

    beam.push_front(bp);
}

// Do not use for the first site! So, callLocIndex has to be >= 1!
// if referenceIsTwoHetsBack is set to 'true', callLocIndex has to be >= 2!
fp_type Beam::callKeepPhaseProb(size_t callLocIndex, bool referenceIsTwoHetsBack) {

    // check bounds (should never be violated)
    if (callLocIndex >= callLocs.size()) {
        StatusFile::addError("callKeepPhaseProb() called with callLocIndex out of bounds!");
        exit(EXIT_FAILURE);
    }

    int frontloc = reverseBeam ?
            max((int)callLocs[callLocIndex] - (int)DeltaCall, 0) :
            min(callLocs[callLocIndex] + DeltaCall, constraints.size()-1); // beam front (actual call loc + Delta) (constraints is as large as splits)
    // where to find the actual haps and the reference haps in the path histories
    curridx = reverseBeam ? callLocs[callLocIndex] - frontloc : frontloc - callLocs[callLocIndex];
    refidx = reverseBeam ?
            (referenceIsTwoHetsBack ? callLocs[callLocIndex+2] : callLocs[callLocIndex+1]) :
            (referenceIsTwoHetsBack ? callLocs[callLocIndex-2] : callLocs[callLocIndex-1]);
    refidx = reverseBeam ? refidx - frontloc : frontloc - refidx;

    // forward beam front to desired position + call length (Delta)
    while((reverseBeam && cursor > frontloc) || (!reverseBeam && cursor < frontloc)) {
        advanceBeamFront();
    }

    // calculate keep against switch probabilities over all actual paths
    // 'keep' and 'switch' is decided according to the actual haplotype pair versus the reference pair

#ifdef DEBUG_TARGET
    cout << " cursor: " << cursor << " curridx: " << curridx << " refidx: " << refidx << endl;
#endif

    // accumulate probabilities
//#ifdef STOPWATCH
//    Stopwatch swcall("callKeepPhase_accProb");
//#endif
//    unsigned count = 0;
    fp_type acckp = 0.0; // keep probability
    fp_type accsp = 0.0; // switch probability
#ifdef DEBUG_TARGET
    int beamsize = 0;
    int keeps = 0;
    int switches = 0;
    int homs = 0;
#endif
    for (auto &bp : beam) {

        auto currhaps = bp.getHapsAt(curridx);
        auto refhaps = bp.getHapsAt(refidx);

        // calculate switch probability against keep probability, whereby the actual haps are compared to their reference one (or two) call locs back
        fp_type p = bp.getPathWeight() * bp.getPathProbability(); // actual weight * actual probability (at beam front)
        if (currhaps.first != currhaps.second && refhaps.first != refhaps.second) {
            // both call locs in the path are heterozygous, so count keep probability against switch probability
            if (currhaps.first == refhaps.first) { // keep
                acckp += p;
#ifdef DEBUG_TARGET
                keeps++;
#endif
            } else { // switch
                accsp += p;
#ifdef DEBUG_TARGET
                switches++;
#endif
            }
        } else { // although callLocs are always heterozygous at the target genotype, the paths might be homozygous at the call locations
            // one of the call locs is homozygous, so we cannot speak of a keep or switch here:
            // according to Eagle, the probability is equally divided to switch and keep probability
            acckp += p/2;
            accsp += p/2;
#ifdef DEBUG_TARGET
            homs++;
#endif
        }
#ifdef DEBUG_TARGET
        beamsize++;
#endif

    }

    fp_type ret = 0.5; // default: fifty-fifty
    // keep probability against switch probability (is 1 if there are no switch paths, regardless of path probability)
    if (acckp + accsp > 0)
        ret = acckp / (acckp + accsp);
#ifdef DEBUG_TARGET
    cout << " beamsize: " << beamsize << " keeps: " << keeps << " " << acckp << " switches: " << switches << " " << accsp << " (homs: " << homs << ")" << endl;
#endif

//#ifdef STOPWATCH
//    swcall.stop();
//#endif
    return ret;
}

fp_type Beam::callDosage(size_t splitLocIndex) {

    // check bounds (should never be violated)
    if (splitLocIndex >= constraints.size()) {
        StatusFile::addError("callDosage() called with splitLocIndex out of bounds!");
        exit(EXIT_FAILURE);
    }

    int frontloc = reverseBeam ?
            max((int)splitLocIndex - (int)DeltaCall, 0) :
            min(splitLocIndex + DeltaCall, constraints.size()-1); // beam front (actual call loc + Delta) (constraints is as large as splits)
    // where to find the actual haps in the path histories
    curridx = reverseBeam ? splitLocIndex - frontloc : frontloc - splitLocIndex;
    refidx = -1; // not required here, set to nothing just for clean debug output

    // forward beam front to desired position
    while((reverseBeam && cursor > frontloc) || (!reverseBeam && cursor < frontloc)) {
        advanceBeamFront();
    }

//#ifdef STOPWATCH
//    Stopwatch swdose("dosage");
//#endif

    fp_type dos = 0.0;
    fp_type ptot = 0.0;
#ifdef DEBUG_TARGET
    int beamsize = 0;
    int homrefs = 0;
    int hets = 0;
    int homalts = 0;
#endif

    // calculate dosage (probability modified by current haps) against total path probability at beam front
    for (auto &bp : beam) {
        auto currhaps = bp.getHapsAt(curridx);

        fp_type p = bp.getPathWeight() * bp.getPathProbability(); // actual weight * actual probability (at beam front)
        ptot += p;
        if (currhaps.first == currhaps.second) { // homozygous
            if (currhaps.first == Haplotype::Alt) { // hom alt -> *2
                dos += 2*p;
#ifdef DEBUG_TARGET
                homalts++;
#endif
            } else { // hom ref -> *0
#ifdef DEBUG_TARGET
                homrefs++;
#endif
            }
        } else { // heterozygous -> *1
            dos += p;
#ifdef DEBUG_TARGET
            hets++;
#endif
        }
#ifdef DEBUG_TARGET
        beamsize++;
#endif
    }

    fp_type ret = 1.0; // default: corresponds to fifty-fifty (diploid dosage!)
    if (ptot > 0.0)
        ret = dos / ptot;

#ifdef DEBUG_TARGET
    cout << "beamsize: " << beamsize << " cursor: " << cursor << " curridx: " << curridx << " homrefs: " << homrefs << " hets: " << hets << " homalts: " << homalts << " (dos: " << dos << " / ptotal: " << ptot << " = " << ret << ")" << endl;
#endif

//#ifdef STOPWATCH
//    swdose.stop();
//#endif

    return ret;
}

void Beam::advanceBeamFront() {
    // ATTENTION! No boundary checks for speed reasons! Be sure not to advance past the end!

//#ifdef STOPWATCH
//    Stopwatch swadvance("advanceBeamFront");
//#endif

    // advance cursor (runs over split sites!)
    if (reverseBeam)
        cursor--;
    else
        cursor++;

    pbwt.advanceTo(2*cursor+1); // the odd sites mark the reference sites, evens are incon segments

    int32_t bpscaleexp = 0; // global beampath scaling exponent (determined later during happath extension)

#ifdef DEBUG_TARGET
    {
        int cnt0 = pbwt.getCount0(2*cursor+1);
        cout << "advance(): cursor: " << cursor
             << " count0: " << cnt0 << " (" << (100*cnt0*pbwt.getKInv()) << "%)";
        cnt0 = pbwt.getCount0(2*(cursor + reverseBeam ? 1 : 0));
        cout << " incons: " << cnt0 << " (" << (100*cnt0*pbwt.getKInv()) << "%)" << endl;
    }
#endif

    // calculate recombination probabilities
//#ifdef STOPWATCH
//    Stopwatch swcr("calcRec");
//#endif
    calculateRecProbs();
//#ifdef STOPWATCH
//    swcr.stop();
//#endif

    bool repeat;
    Constraint constr = constraints[cursor];
    do {
        // 1. Extend current beam paths and
        // 2. compute new path probabilities
//    #ifdef STOPWATCH
//        Stopwatch swex("extendPaths");
//    #endif
        fp_type maxP = 0.0; // maximum probability in beam
        {
    #ifdef DEBUG_MAPPINGS_TO_STDERR
            vector<DbgMapContainer> dbgmappings;
    #endif
//    #ifdef STOPWATCH
//            Stopwatch swmark("mark ext");
//    #endif
            // firstly, iterate through beam and mark all required happath extensions
            for (auto &bp : beam) {
                bp.markRequiredExtensions(constr);
            }
//    #ifdef STOPWATCH
//            swmark.stop();
//    #endif

//    #ifdef STOPWATCH
//            Stopwatch swhp("hap prob");
//    #endif
            // process extensions in happaths and calculate happath probabilities
            const HaplotypePath *hp_neighb = nullptr; // null
            int32_t minhpscaleexp = 0;
            int32_t maxhpscaleexp = INT32_MIN;
            for (size_t hpi = 0; hpi < happaths.size(); hpi++) {
                HaplotypePath *hp = happaths[hpi];
    #if defined DEBUG_TARGET && defined DEBUG_PRINT_PBWT_MAPPINGS
                cout << "happath: " << hpi << " " << hp->getPathProbability().p << " " << hp->getPathProbability().scaleexp << endl;
    #endif
                if (hp->isToBeExtended()) {
                    int32_t scaleexp = hp->extend(prec, cursor, hp_neighb
    #if defined DEBUG_TARGET && defined DEBUG_PRINT_PBWT_MAPPINGS
                     , constr
    #endif
    #ifdef DEBUG_MAPPINGS_TO_STDERR
                     , dbgmappings
    #endif
                     );
                    hp_neighb = hp; // save current path as neighbor for the next extension
                    if (scaleexp < minhpscaleexp)
                        minhpscaleexp = scaleexp;
                    if (scaleexp > maxhpscaleexp)
                        maxhpscaleexp = scaleexp;
    #if defined DEBUG_TARGET && defined DEBUG_PRINT_PBWT_MAPPINGS
                    cout << "  p0: " << hp->getExtensionProbability(0).p << " " << hp->getExtensionProbability(0).scaleexp
                         << (hp->isExtensionActive(0) ? " active" : " inactive")
                         << "  p1: " << hp->getExtensionProbability(1).p << " " << hp->getExtensionProbability(1).scaleexp
                         << (hp->isExtensionActive(1) ? " active" : " inactive") << endl;
    #endif
                }
            }
            // this scaling normalizes the product of the most and the least probable happaths
            bpscaleexp = -(maxhpscaleexp + minhpscaleexp);
//    #ifdef STOPWATCH
//            swhp.stop();
//    #endif

//    #ifdef STOPWATCH
//            Stopwatch swpp("path prob");
//    #endif

            // calculate beam path probabilities
            for (auto &bp : beam) {
                fp_type maxPtmp = bp.extend(target[cursor], constr, bpscaleexp);

                if (maxP < maxPtmp)
                    maxP = maxPtmp;
            }
//    #ifdef STOPWATCH
//            swpp.stop();
//    #endif
    #ifdef DEBUG_MAPPINGS_TO_STDERR
            // sort mappings
            sort(dbgmappings.begin(), dbgmappings.end());
            // output only unique entries, but check if "equal" entries are really equal!
            cerr << "cursor: " << cursor+1 << endl;
            auto dbgit = dbgmappings.cbegin();
            while (dbgit != dbgmappings.cend()) {
                for (int h = min(dbgit->history_length-1,63); h >= 0; h--) {
                    unsigned long long mask = 1ull << h;
                    if (dbgit->history & mask)
                        cerr << "1";
                    else
                        cerr << "0";
                }
                cerr << ": " << dbgit->mapping[0] << " " << dbgit->mapping[1] << endl;
                const auto previt = dbgit;
                dbgit++;
                bool goon = true;
                while (dbgit != dbgmappings.end() && *dbgit == *previt && goon) {
                    // equal neighbors -> check if they are really equal!
                    if (dbgit->mapping[0] == previt->mapping[0] && dbgit->mapping[1] == previt->mapping[1])
                        // yes, go on!
                        dbgit++;
                    else {
                        // no! Oh oh, this should never be the case!
                        cerr << "!!! ";
                        goon = false; // break
                    }
                }
            }
            cerr << endl;
    #endif
        }
//    #ifdef STOPWATCH
//        swex.stop();
//    #endif

    #if defined DEBUG_TARGET && defined DEBUG_PRINT_PATHS_AFTER_EXT
        {
            cout << "paths after extension: scaleexp: " << bpscaleexp << endl;
            int bpi = 0;
    #ifdef DEBUG_PRINT_SORTED_BY_PROB
            list<BeamPath> beamcpy((const list<BeamPath>)beam); // copy
            beamcpy.sort();
    #else
            list<BeamPath> &beamcpy = beam;
    #endif
            for (auto &bp : beamcpy) {
                cout << " " << bpi << ": " << bp.getPathWeight() << " " << bp.getPathProbability() << endl;
                cout << "  ";
                bp.getPathHistoryMat().printPath(cout, curridx-1);
                cout << "  " << bp.getMatPathProb().p << " " << bp.getMatPathProb().scaleexp;
                cout << "\n  ";
                bp.getPathHistoryPat().printPath(cout, curridx-1);
                cout << "  " << bp.getPatPathProb().p << " " << bp.getPatPathProb().scaleexp;
                cout << endl;
                for (int j = 0; j < 4; j++) {
                    cout << "  ext " << j << ": " << bp.getExtensionProbability(j)
                         << "\t(m: ";
                    if (bp.isExtMatPathActive(j))
                        cout << bp.getExtMatPathProb(j).p << " " << bp.getExtMatPathProb(j).scaleexp;
                    else
                        cout << "---";
                    cout << " p: ";
                    if (bp.isExtPatPathActive(j))
                        cout << bp.getExtPatPathProb(j).p << " " << bp.getExtPatPathProb(j).scaleexp;
                    else
                        cout << "---";
                    cout << ")\t" << (bp.isExtensionActive(j) ? " active" : " inactive") << endl;
                }
                bpi++;
            }
        }
    #endif

//    #ifdef STOPWATCH
//        Stopwatch swprune("prunePaths");
//    #endif
        // 3. prune probabilities to maintain beam width
        if (prunePaths(maxP, constr))
            repeat = false; // prune successful, beam is not empty
        else { // prune not successful -> empty beam!
            // if there was a constraint applied, we remove the constraint and try again
            if (constr.type != Constraint::Type::NoConstr) {
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
                cout << "Beam is empty! Trying again without constraint." << endl;
#endif
                constr.type = Constraint::Type::NoConstr;
                repeat = true;
            } else { // nothing we could do... TODO better way to exit than crashing?
                cout << endl;
                StatusFile::addError("Beam is empty. No paths left.");
                exit(EXIT_FAILURE);
            }
        }
//    #ifdef STOPWATCH
//        swprune.stop();
//    #endif
    } while(repeat);

#if defined DEBUG_TARGET && defined DEBUG_PRINT_PATHS_AFTER_EXT_APP
    {
        cout << "paths after extension application: scaleexp: " << bpscaleexp << endl;
        int bpi = 0;
#ifdef DEBUG_PRINT_SORTED_BY_PROB
        list<BeamPath> beamcpy((const list<BeamPath>)beam); // copy
        beamcpy.sort();
#else
        list<BeamPath> &beamcpy = beam;
#endif
        for (auto bpit = beamcpy.begin(); bpit != beamcpy.end(); bpit++, bpi++) {
            cout << " " << bpi << ": " << bpit->getPathWeight() << " " << bpit->getPathProbability() <<
                    " #err: " << bpit->getPathNumErrors() << " (#forgot: " << bpit->getPathNumForgottenErrors() << ")" <<
                    " pErr: " << bpit->getPathErrorProb() << endl;
            cout << "  ";
            bpit->getPathHistoryMat().printPath(cout, curridx, refidx);
            cout << "  " << bpit->getMatPathProb().p << " " << bpit->getMatPathProb().scaleexp;
            cout << "\n  ";
            bpit->getPathHistoryPat().printPath(cout, curridx, refidx);
            cout << "  " << bpit->getPatPathProb().p << " " << bpit->getPatPathProb().scaleexp;
            cout << endl;
//            cout << "  p/m: " << (bpit->getMatPathProb() / bpit->getPatPathProb()) << endl;
        }
    }
#endif

//#ifdef STOPWATCH
//    Stopwatch swmerge("mergePaths");
//#endif
    // 4. Merge paths
    mergePaths();
//#ifdef STOPWATCH
//    swmerge.stop();
//#endif

#if defined DEBUG_TARGET && defined DEBUG_PRINT_PATHS_AFTER_MERGE
    {
        cout << "paths after merge: scaleexp: " << bpscaleexp << endl;
        int bpi = 0;
#ifdef DEBUG_PRINT_SORTED_BY_PROB
        list<BeamPath> beamcpy((const list<BeamPath>)beam); // copy
        beamcpy.sort();
#else
        list<BeamPath> &beamcpy = beam;
#endif
        for (auto &bp : beamcpy) {
            cout << " " << bpi << ": " << bp.getPathWeight() << " " << bp.getPathProbability() << " hpi: " << bp.getPathMatIndex() << " + " << bp.getPathPatIndex() << endl;
            cout << "  ";
            bp.getPathHistoryMat().printPath(cout, curridx-1);
            cout << "  " << bp.getMatPathProb().p << " " << bp.getMatPathProb().scaleexp;
            cout << "\n  ";
            bp.getPathHistoryPat().printPath(cout, curridx-1);
            cout << "  " << bp.getPatPathProb().p << " " << bp.getPatPathProb().scaleexp;
            cout << endl;
            bpi++;
        }
    }
#endif

//#ifdef STOPWATCH
//    swadvance.stop();
//#endif
}

void Beam::calculateRecProbs() {

    int c = reverseBeam ? cursor + 1 : cursor; // running cursor for cMsplitDiffs (to the left)
    fp_type cursordiff = reverseBeam ? cMsplitDiffs[cursor] : cMsplitDiffs[cursor+1];

    fp_type curr_diff = 0.0; // this is the case if split and cursor are at the same position (beginning of loop)
    for (int prec_i = 0; prec_i < (int)H; prec_i++) { // most recent recombination probability is always at index 0
        // default in Eagle
        fp_type tmp1 = 1.0 + curr_diff / expectIBDcM;
        fp_type pr = 1.0 / (tmp1 * tmp1);
        // cMsplitDiffs is padded, so the info can be used also on last site (forward) or first site (reverse)
        fp_type tmp2 = 1.0 + (curr_diff + cursordiff) / expectIBDcM; // use difference of current cursor site to next site to get d_{x:m+1}
        pr -= 1.0 / (tmp2 * tmp2);
        pr = max(min(pr, (fp_type)1.0), (fp_type)0.000001); // ensure 0 < pr <= 1

#if defined DEBUG_TARGET && defined DEBUG_PRINT_PREC
        cout << " prec " << prec_i << ": " << pr << endl;
#endif
        prec[prec_i] = pr * pbwt.getKInv(); // P^(rec m|x) = (1/Kc_x) * P(rec m|x)

        curr_diff += cMsplitDiffs[c]; // add difference from last to current position to get d_{c:cursor}

        if (reverseBeam) {
            c++; // running cursor to the right for reverse Beam
            if (c >= (int)cMsplitDiffs.size())
                break;
        } else {
            if (c == 0)
                break;
            c--; // running cursor to the left
        }
    }
}



bool Beam::prunePaths(fp_type maxP, const Constraint &constr) {
//#ifdef STOPWATCH
//    Stopwatch swfilter("filter");
//#endif

    // filter and sort according to probability
    fp_type limit = pLimit * maxP;
    MinMaxHeap<BeamPathContainer> beamProbs;

#if defined DEBUG_TARGET && defined DEBUG_PRINT_PATHS_AFTER_PRUNE
    cout << "maxP: " << maxP << " prune limit: " << limit << endl;
    fp_type totP = 0.0, totExP = 0.0;
#endif
    for (beam_iterator it = beam.begin(); it != beam.end(); it++) {
#if defined DEBUG_TARGET && defined DEBUG_PRINT_PATHS_AFTER_PRUNE
        totP += it->getPathProbability() * it->getPathWeight();
#endif
        int ex = it->getActiveExtensionIndex();
        while (ex >= 0) { // over all possible extension indices
            fp_type p = it->getExtensionProbability(ex);
#if defined DEBUG_TARGET && defined DEBUG_PRINT_PATHS_AFTER_PRUNE
//                    totExP += p;
            totExP += it->getPathWeight() * p;
#endif
            if (p >= limit) { // element will potentially be kept
                if (beamProbs.getSize() < P || (beamProbs.getSize() == P && p > beamProbs.getMin().prob)) {
                    // insert into MinMaxHeap
                    if (beamProbs.getSize() == P) {
                        // remove smallest element
                        BeamPathContainer bpc = beamProbs.popMin();
                        bpc.it->markExtensionDeleted(bpc.extensionIndex);
                    }
                    beamProbs.emplace(it, p, ex);
                } else // not kept -> mark as deleted
                    it->markExtensionDeleted(ex);
            } else // p < limit -> mark as deleted
                it->markExtensionDeleted(ex);
            // only continue if we are unconstrained (only then all extensions are active, otherwise it's only 1)
            ex = constr.type == Constraint::Type::NoConstr ? ex-1 : -1;
        } // END over all extension indices
    } // END for all paths

#if defined DEBUG_TARGET && defined DEBUG_PRINT_PATHS_AFTER_PRUNE
    {
        cout << "paths after prune:" << endl;
        int bpi = 0;
#ifdef DEBUG_PRINT_SORTED_BY_PROB
        list<BeamPath> beamcpy((const list<BeamPath>)beam); // copy
        beamcpy.sort();
#else
        list<BeamPath> &beamcpy = beam;
#endif
        for (auto &bp : beamcpy) {
            fp_type p = bp.getPathWeight() * bp.getPathProbability();
            cout << " " << bpi << ": " << bp.getPathWeight() << " x " << bp.getPathProbability() << " = " << p << " (" << p/totP << ")" << endl;
            cout << "  ";
            bp.getPathHistoryMat().printPath(cout, curridx-1);
            cout << "  " << bp.getMatPathProb().p << " " << bp.getMatPathProb().scaleexp;
            cout << "\n  ";
            bp.getPathHistoryPat().printPath(cout, curridx-1);
            cout << "  " << bp.getPatPathProb().p << " " << bp.getPatPathProb().scaleexp;
            cout << endl;
//            cout << "  p/m: " << (bp.getMatPathProb() / bp.getPatPathProb()) << endl;
            for (int j = 0; j < 4; j++)
                cout << "  ext " << j << ": " << bp.getExtensionProbability(j) << "\t(" << (bp.getPathWeight() * bp.getExtensionProbability(j))/totExP << ")"
                 << "\t(m: " << bp.getExtMatPathProb(j).p << " " << bp.getExtMatPathProb(j).scaleexp
                 << " p: " << bp.getExtPatPathProb(j).p << " " << bp.getExtPatPathProb(j).scaleexp << ")"
                 << "\t" << (bp.isExtensionActive(j) ? " active" : (bp.getExtensionProbability(j) >= limit ? " inactive (by order)" : " inactive (by limit/init)")) << endl;
//            cout << "  m-ext 0/1: ";
//            if (bp.getExtMatPathProb(0) > 0.0 && bp.getExtMatPathProb(2) > 0.0)
//                cout << (bp.getExtMatPathProb(0) / bp.getExtMatPathProb(2));
//            else
//                cout << "na";
//            cout << "  p-ext 0/1: ";
//            if (bp.getExtPatPathProb(0) > 0.0 && bp.getExtPatPathProb(1) > 0.0)
//                cout << (bp.getExtPatPathProb(0) / bp.getExtPatPathProb(1));
//            else
//                cout << "na";
//            cout << endl;
            bpi++;
        }
    }
#endif

//#ifdef STOPWATCH
//    swfilter.stop();
//#endif

    // now, the MinMaxHeap has kept track of all elements that should remain
    // and all elements, that should be deleted are marked.
    // so, iterate through the beam and apply the active extensions as well as erase marked elements

    // if the heap is empty, return false as an unsuccessful prune
    if (beamProbs.isEmpty())
        return false;

//#ifdef STOPWATCH
//    Stopwatch swapply("apply");
//#endif

    // we could also apply a simple weight scaling here: we take the exponent of the weight of the beam with the smallest probability as (negative) scaling factor
    const BeamPathContainer &minpath = beamProbs.getMin();
    int weightScaleExp = -ilogb(minpath.it->getPathWeight());
    // and we do error scaling by letting all paths forget an error if the most probable path contains at least one
    // (we assume that the most probable path always contains the smallest number of errors, but anyway, it won't be a problem if not)
    bool deleteErr = beamProbs.getMax().it->getPathNumErrors() > 0;

//#ifdef STOPWATCH
//    Stopwatch swapplyhap("apply hap");
//#endif

    // apply haplotype path extensions
    vector<HaplotypePath*> tmp_happaths = happaths; // copy
    int nhapext0 = 0;
    for (const auto hap : happaths) {
        // count 0-extensions, but not twice for double entries
        if (hap->isExtensionActive(0))// && hap != hp_neighb)
            nhapext0++;
    }
    vector<pair<int,int>> hapmap; // map which index maps where
    hapmap.reserve(happaths.size());
    happaths.resize(nhapext0); // this is the minimum new size (if there were no 1-extensions)
    int curroneidx = nhapext0; // the index of the first 1-extension
    int currzeroidx = 0;
    for (const auto tmphap : tmp_happaths) { // over all old happaths

        // mapping: if one path is not extended the following is still ok since the corresponding value will not be used.
        hapmap.push_back(make_pair(currzeroidx, curroneidx));

        // no active extension -> delete object
        if (tmphap->allExtensionsInactive()) {
            // object is not required anymore
            delete tmphap;
            continue;
        }

        // all extensions active
        if (tmphap->allExtensionsActive()) {
            // copy path for one-extension, apply zero-extension
            happaths.push_back(new HaplotypePath(*tmphap, 1));
            curroneidx++;
            tmphap->applyExtension(0);
            happaths[currzeroidx] = tmphap;
            currzeroidx++;
            continue;
        }

        // here either one- or zero-extension is active -> apply it
        // apply zero-extension (and copy element to new position)
        if (tmphap->isExtensionActive(0)) {
            tmphap->applyExtension(0);
            happaths[currzeroidx] = tmphap;
            currzeroidx++;
        } else {
            tmphap->applyExtension(1);
            happaths.push_back(tmphap);
            curroneidx++;
        }
    }

//#ifdef STOPWATCH
//    swapplyhap.stop();
//#endif
//
//#ifdef STOPWATCH
//    Stopwatch swapplybeam("apply beam");
//#endif

    // save the iterator to the last element
    beam_iterator lastit = beam.end();
    lastit--;

    beam_iterator bpit = beam.begin();
    do {
        int ex = bpit->getActiveExtensionIndex();
        if (ex >= 0) { // there's at least one active extension
            // getting the mapping indices
            int mat[2], pat[2];
            mat[0] = hapmap[bpit->getPathMatIndex()].first;
            mat[1] = hapmap[bpit->getPathMatIndex()].second;
            pat[0] = hapmap[bpit->getPathPatIndex()].first;
            pat[1] = hapmap[bpit->getPathPatIndex()].second;

            // for all active extensions but the last (that will be applied) generate a new BeamPath object
            int next_ex = bpit->getNextActiveExtensionIndex(ex);
            while (next_ex != -1) { // there will still be another active extension
                beam.emplace_back(*bpit, mat[ex/2], pat[ex%2], ex, weightScaleExp, deleteErr); // insert element at the end of the list
                ex = next_ex;
                next_ex = bpit->getNextActiveExtensionIndex(next_ex);
            }
            // the last active extension will be applied
            bpit->applyExtension(mat[ex/2], pat[ex%2], ex, weightScaleExp, deleteErr);
            if (bpit == lastit)
                break;
            else
                bpit++;
        } else { // no active extensions -> delete element
            if (bpit == lastit) {
                beam.erase(bpit);
                break;
            } else
                bpit = beam.erase(bpit); // now bpit points to next element
        }
    } while(1);

//#ifdef STOPWATCH
//    swapplybeam.stop();
//    swapply.stop();
//#endif
    return true;
}

// merge paths that have equal haplotype history at the last history sites before extension!
void Beam::mergePaths() {

    unordered_map<size_t, beam_iterator> mergemap;
    beam_iterator bpit = beam.begin();
    while (bpit != beam.end()) {
        size_t matref = bpit->getPathMatIndex();
        // find equal happath with smaller index
        while (matref > 0 && happaths[matref]->isMergeable(*(happaths[matref-1])))
            matref--;
        size_t patref = bpit->getPathPatIndex();
        // same for pat
        while (patref > 0 && happaths[patref]->isMergeable(*(happaths[patref-1])))
            patref--;
        // make key for map
        size_t k0 = matref;
        size_t k1 = patref;
        if (k0 > k1) { // ensure k0 is smaller or equal than k1
            k0 = patref;
            k1 = matref;
        }
        size_t key = (k0 << 32) | k1; // maybe other compositions for a key work better?
        // look if there's already an element with the same key in the map
        auto mit = mergemap.find(key);
        if (mit != mergemap.end()) { // found an element -> merge
            beam_iterator bpmit = (*mit).second;
            if (bpmit->getPathProbability() > bpit->getPathProbability()) { // element in map is better
#if defined DEBUG_TARGET && defined DEBUG_PRINT_PATHS_AFTER_MERGE
                cout << "Merge: hps " << bpit->getPathMatIndex() << " + " << bpit->getPathPatIndex() << " -> " << bpmit->getPathMatIndex() << " + " << bpmit->getPathPatIndex() << endl;
#endif
                bpmit->merge(*bpit);
                bpit = beam.erase(bpit); // erase element, returns iterator to next element
                // can continue with next element here
            } else { // new element is better
#if defined DEBUG_TARGET && defined DEBUG_PRINT_PATHS_AFTER_MERGE
                cout << "Merge: hps " << bpmit->getPathMatIndex() << " + " << bpmit->getPathPatIndex() << " -> " << bpit->getPathMatIndex() << " + " << bpit->getPathPatIndex() << endl;
#endif
                bpit->merge(*bpmit);
                beam.erase(bpmit);   // delete old path from beam
                mergemap.erase(mit); // delete old path from map
                mergemap[key] = bpit; // insert new one in map
                bpit++; // next element and continue
            }
        } else { // not in map -> insert
            mergemap[key] = bpit;
            bpit++;
        }
    }
}
