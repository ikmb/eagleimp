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

#include <cmath>
#include <algorithm>

#include "Datatypes.h"

#include "BeamPath.h"

BeamPath::BeamPath(fp_type pErr_, vector<HaplotypePath*> &happaths_, int mat, int pat, fp_type weight_, fp_type pErrPath_, int nErrPath_)
 : happaths(happaths_), pathmat(mat), pathpat(pat),
   weight(weight_),
   pErr(pErr_), nErrPath(nErrPath_), pErrPath(pErrPath_), pErrPathp1(pErrPath_*pErr_),
#ifdef DEBUG_TARGET
   forgottenerr(0),
#endif
   activeExtensions{false,false,false,false},
   activeExtensionIndex(-1),
//   activeExtensionCount(0),
   extTargetGeno(Genotype::Miss) // default
{
    HaplotypePathProbability pm = happaths[mat]->getPathProbability();
    HaplotypePathProbability pp = happaths[pat]->getPathProbability();
    prob = pow(2, pm.scaleexp + pp.scaleexp) * pm.p * pp.p * pErrPath_;
    // clear extension fields for new four possible extensions
    extensionProbs.fill(0.0);
}

BeamPath::BeamPath(const BeamPath &other, int mat, int pat, int extensionIndex, int weightScaleExp, bool deleteErr)
 : happaths(other.happaths), pathmat(mat), pathpat(pat),
   prob(other.extensionProbs[extensionIndex]), // set the probability to that of the selected extension
   weight(other.weight), // set the weight of the other path
   pErr(other.pErr), nErrPath(other.nErrPath), pErrPath(other.pErrPath), pErrPathp1(other.pErrPathp1),
#ifdef DEBUG_TARGET
   forgottenerr(other.forgottenerr),
#endif
   activeExtensions{false,false,false,false},
   activeExtensionIndex(-1),
//   activeExtensionCount(0),
   extTargetGeno(other.extTargetGeno)
{
    // apply weight extension, if desired
    if (weightScaleExp != 0)
        weight *= pow(2.0, weightScaleExp);

    // keep track of the current number of errors in path
    updateErrors(extensionIndex, deleteErr);

    // clear extension fields for new four possible extensions
    extensionProbs.fill(0.0);

    extTargetGeno = Genotype::Miss;

}

void BeamPath::applyExtension(int mat, int pat, int extensionIndex, int weightScaleExp, bool deleteErr) {

    pathmat = mat;
    pathpat = pat;

    // set the probability to that of the selected extension
    prob = extensionProbs[extensionIndex];

    // apply weight extension, if desired
    if (weightScaleExp != 0)
        weight *= pow(2.0, weightScaleExp);

    // keep track of the current number of errors in path
    updateErrors(extensionIndex, deleteErr);

    // clear extension fields for new four possible extensions
    extensionProbs.fill(0.0);

    activeExtensions[0] = false;
    activeExtensions[1] = false;
    activeExtensions[2] = false;
    activeExtensions[3] = false;
    activeExtensionIndex = -1;
//    activeExtensionCount = 0;
    extTargetGeno = Genotype::Miss;
}

void BeamPath::merge(BeamPath &other) {
    weight += other.weight * (other.prob / prob);
}

void BeamPath::markExtensionDeleted(int extensionIndex) {
//    if (activeExtensions[extensionIndex]) {
        activeExtensions[extensionIndex] = false;
//        activeExtensionCount--; // assumed to be called only on an active extension!

        // remove the reference from this beam path to the corresponding happaths
        happaths[pathmat]->markExtensionDeleted(extensionIndex/2);
        happaths[pathpat]->markExtensionDeleted(extensionIndex%2);

        updateActiveExtensionIndex();
//    }
}

void BeamPath::markRequiredExtensions(const Constraint &constr) {
    happaths[pathmat]->markRequiredExtensions(constr);
    happaths[pathpat]->markRequiredExtensions(constr);
}

//// returns the weighted maximum probability of all extensions
// returns the maximum probability of all extensions (unweighted)
fp_type BeamPath::extend(Genotype targetGeno, const Constraint &constr, int32_t bpscaleexp) {
    // calculate all probabilities for the four possible extensions with the probabilities of the already extended haplotype paths

    extTargetGeno = targetGeno; // save the genotype in the target at the extension site

    HaplotypePath *pm = happaths[pathmat];
    HaplotypePathProbability pm0 = pm->getExtensionProbability(0);
    HaplotypePathProbability pm1 = pm->getExtensionProbability(1);
    bool calcpm0 = pm->isExtensionActive(0);
    bool calcpm1 = pm->isExtensionActive(1);

    HaplotypePath *pp = happaths[pathpat];
    HaplotypePathProbability pp0 = pp->getExtensionProbability(0);
    HaplotypePathProbability pp1 = pp->getExtensionProbability(1);
    bool calcpp0 = pp->isExtensionActive(0);
    bool calcpp1 = pp->isExtensionActive(1);

    // required scaling for these two paths, adjusted by the global beampath scaling factor
    // (scaling factor is the same for pm0 and pm1, for pp0 and pp1 respectively)
    fp_type hpscaleexp = pow(2, pm0.scaleexp + pp0.scaleexp + bpscaleexp);

    // total path probability for the four possible extensions
    // only calculate those required by the constraint and those with a probability > 0.0
    // and determine maximum
    fp_type maxP = 0.0;
    bool incon = inconsistent(0, targetGeno); // check if the target genotype is inconsistent with this extension
    if (constr.type == Constraint::Type::NoConstr || constr.type == Constraint::Type::HomRef) {
        // extensions for 00 (pm0*pp0)
        if (calcpm0 && calcpp0) {
            extendSingle(0, pm0.p, pp0.p, hpscaleexp, incon, maxP);
            happaths[pathmat]->incExtReferenceCnt(0);
            happaths[pathpat]->incExtReferenceCnt(0);
        }
    }
    incon = inconsistent(1, targetGeno); // check if the target genotype is inconsistent with this extension
    if (constr.type == Constraint::Type::NoConstr
            || (constr.type == Constraint::Type::HetKeep && pm->getPathHistory()[constr.refidx] == Haplotype::Ref)
            || (constr.type == Constraint::Type::HetSwitch && pm->getPathHistory()[constr.refidx] == Haplotype::Alt)) {
        // extensions for 01 (pm0*pp1)
        if (calcpm0 && calcpp1) {
            extendSingle(1, pm0.p, pp1.p, hpscaleexp, incon, maxP);
            happaths[pathmat]->incExtReferenceCnt(0);
            happaths[pathpat]->incExtReferenceCnt(1);
        }
    }
    incon = inconsistent(2, targetGeno); // check if the target genotype is inconsistent with this extension
    if (constr.type == Constraint::Type::NoConstr
            || (constr.type == Constraint::Type::HetKeep && pm->getPathHistory()[constr.refidx] == Haplotype::Alt)
            || (constr.type == Constraint::Type::HetSwitch && pm->getPathHistory()[constr.refidx] == Haplotype::Ref)) {
        // extensions for 10 (pm1*pp0)
        if (calcpm1 && calcpp0) {
            extendSingle(2, pm1.p, pp0.p, hpscaleexp, incon, maxP);
            happaths[pathmat]->incExtReferenceCnt(1);
            happaths[pathpat]->incExtReferenceCnt(0);
        }
    }
    incon = inconsistent(3, targetGeno); // check if the target genotype is inconsistent with this extension
    if (constr.type == Constraint::Type::NoConstr || constr.type == Constraint::Type::HomAlt) {
        // extensions for 11 (pm1*pp1)
        if (calcpm1 && calcpp1) {
            extendSingle(3, pm1.p, pp1.p, hpscaleexp, incon, maxP);
            happaths[pathmat]->incExtReferenceCnt(1);
            happaths[pathpat]->incExtReferenceCnt(1);
        }
    }

//    return weight * maxP;
    return maxP;
}

inline void BeamPath::extendSingle(int id, fp_type p1, fp_type p2, fp_type scale, bool incon, fp_type &maxP) {
    extensionProbs[id] = p1 * p2 * scale;
    if (incon)
        extensionProbs[id] *= pErrPathp1;
    else
        extensionProbs[id] *= pErrPath;
    if (maxP < extensionProbs[id])
        maxP = extensionProbs[id];
    activeExtensions[id] = true;
    activeExtensionIndex = id; // assumed that this function is called in ascending order for all extensions
//    activeExtensionCount++; // assumed not to be called twice with the same index
}

inline void BeamPath::updateErrors(int extensionIndex, bool deleteErr) {
    // the number of errors to be added (or removed if negative)
    int adderr = errorsToAdd(extensionIndex, extTargetGeno) - (deleteErr ? 1 : 0);
    nErrPath += adderr;

    if (adderr < 0) { // negative, we have to forget one error (could only be -1)
        pErrPathp1 = pErrPath;
        pErrPath /= pErr;
    } else { // positive or zero
        for (int a = 0; a < adderr; a++) { // decrease the total path error probability by the number of errors to be added
            pErrPath = pErrPathp1;
            pErrPathp1 *= pErr;
        }
    }
#ifdef DEBUG_TARGET
    forgottenerr += (deleteErr ? 1 : 0); // save the number of errors we forgot during analysis
#endif
}

int BeamPath::getNextActiveExtensionIndex(int x) {
    int ret = -1;
    for (int i = x-1; i >= 0; i--) {
        if (activeExtensions[i])
            return i;
    }
    return ret;
}

void BeamPath::updateActiveExtensionIndex() {
    activeExtensionIndex = -1;
    for (int i = 3; i >= 0; i--) {
        if (activeExtensions[i]) {
            activeExtensionIndex = i;
            break;
        }
    }
}

inline pair<Haplotype,Haplotype> BeamPath::getHapsFromExtensionIndex(int extensionIndex) const {
    int extInd = extensionIndex;
    switch (extInd) {
    case 0:
        return make_pair(Haplotype::Ref,Haplotype::Ref);
    case 1:
        return make_pair(Haplotype::Ref,Haplotype::Alt);
    case 2:
        return make_pair(Haplotype::Alt,Haplotype::Ref);
    default: // 3
        return make_pair(Haplotype::Alt,Haplotype::Alt);
    }
}

inline bool BeamPath::inconsistent(int extInd, Genotype targetGeno) const {
    switch (targetGeno) {
    case Genotype::Het:
        if (extInd == 0 || extInd == 3) {
            return true;
        }
        break;
    case Genotype::HomRef:
        if (extInd > 0) {
            return true;
        }
        break;
    case Genotype::HomAlt:
        if (extInd < 3) {
            return true;
        }
        break;
    default: // missing, should never occur here
        return true;
    }
    return false;
}

inline int BeamPath::errorsToAdd(int extensionIndex, Genotype targetGeno) const {
    return inconsistent(extensionIndex, targetGeno) ? 1 : 0;
}

bool BeamPath::operator<(const BeamPath &other) const {
//    return prob * weight > other.prob * other.weight;
    return prob > other.prob;
}

