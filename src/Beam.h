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

#ifndef BEAM_H_
#define BEAM_H_

#include <utility>
#include <vector>
#include <list>

#include "History.h"
#include "Datatypes.h"
#include "PBWT.h"
#include "HaplotypePath.h"
#include "BeamPath.h"

using namespace std;

using beam_iterator = list<BeamPath>::iterator;

class Beam {

public:
    Beam(size_t historySize, size_t beamWidth, size_t deltaCall,
            fp_type pErr, fp_type pLimit, fp_type expectIBDcM,
            const vector<fp_type> &cMsplitDiffs_, const vector<size_t> &callLocs,
            const GenotypeVector &target,
            PBWT &pbwt,
            const vector<Constraint> &constraints,
            bool reverseBeam = false);
    ~Beam();

    // returns probability of keeping phase from call location at callLocs[callLocIndex-1] to callLocs[callLocIndex]
    // (so not applicable for callLocIndex == 0)!
    // if referenceIsTwoHetsBack is set to true, the phase probability is called from callLocs[callLocIndex-2] to callLocs[callLocIndex]
    // (so not applicable for callLocIndex <= 1)!
    // reverse beam: +1 resp. +2 instead of -1 resp. -2
    fp_type callKeepPhaseProb(size_t callLocIndex, bool referenceIsTwoHetsBack = false);
    // call dosage for homozygous sites
    fp_type callDosage(size_t splitLocIndex);

private:

    void initBeamFront();
    void advanceBeamFront();
    void calculateRecProbs();

    // keep only the best P paths with probability greater zero, apply extensions
    // requires current maximum path probability and actual constraint,
    // returns true if successful, i.e. beam is not empty
    bool prunePaths(fp_type maxP, const Constraint & constr);
    // merge paths that have equal haplotype history
    void mergePaths();

    size_t H; // history sites
    size_t P; // beam width
    size_t DeltaCall; // call length
    fp_type pErr; // switch error probability
    fp_type pLimit; // probability deviation limit in beam (usually switch error probability or squared sw err prob)
    fp_type expectIBDcM; // expected IBD length in cM

    /* as in PBWTPhaser.h */
    const vector<fp_type> &cMsplitDiffs; // difference between split sites in cM (cMsplitDiffs[1] = cMpos[splits[1]] - cMpos[splits[0]],...) (cMsplitDiffs[0] = cMpos[splits[0]] - cMpos[0]) (cMsplitDiffs[splits.size()] = cMpos.back() - cMpos[splits.back()])
    const vector<size_t> &callLocs; // indices of call locations (i.e. het + not mono) pointing into _condensed_ structure (e.g. target, reference, splits)
    const GenotypeVector &target;   // condensed genotype information of target only at split sites
    const vector<Constraint> &constraints; // extension constraints for each split site

    int cursor; // current split site index at the beam front
    int curridx;
    int refidx;
    // PBWT structure (reference: object provided in constructor)
    PBWT &pbwt;
    // the beam paths including probs and weights of the beam at the current cursor position,
    // each element keeping two pointers to the corresponding haplotype paths
    // (kept at max. P elements and sorted by their haplotype history)
    list<BeamPath> beam;
    // haplotype paths sorted by their haplotype history
    vector<HaplotypePath*> happaths;
    // recombination probabilities prepared for next extension of BeamPath
    // (independent of BeamPath, so calculated here), size kept at H although not all might be used
    vector<fp_type> prec;
    // true if this beam is evolved in the reverse direction
    bool reverseBeam;

    // if the maximum probability's exponent drops below the limit, we increase the scaling factor by this value
    static const int minProbExp;

};

class BeamPathContainer {
public:
    beam_iterator it;
    fp_type prob;
    int extensionIndex;
    BeamPathContainer(beam_iterator it_, fp_type prob_, int extensionIndex_) : it(it_), prob(prob_), extensionIndex(extensionIndex_) {}
};

#endif /* BEAM_H_ */
