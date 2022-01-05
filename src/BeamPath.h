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

#ifndef BEAMPATH_H_
#define BEAMPATH_H_

#include <utility>
#include <memory>
#include <array>
#include <vector>

#include "History.h"
#include "Datatypes.h"
#include "HaplotypePath.h"
#include "PBWT.h"

class BeamPath {

public:
    // creates a beam path with general error probability pErr, consisting of both provided haplotype paths, weight and current path error probability and number of errors in path
    BeamPath(fp_type pErr, vector<HaplotypePath*> &happaths, int mat, int pat, fp_type weight, fp_type pErrPath, int nErrPath);
    // special constructor that copies the given BeamPath, but applies the selected (and prepared) extension right away
    // the selected extension in the other object will become invalid afterwards
    // also applies an extension scaling, if provided != 0
    // (also includes the error adjustment (scaling), i.e. if 'deleteErr' is true, one error in the path is removed ('forgotten'))
    BeamPath(const BeamPath& other, int mat, int pat, int extensionIndex, int weightScaleExp, bool deleteErr);
    // default copy + move constructors, copy + move assignments
    BeamPath(const BeamPath&) = default;
    BeamPath(BeamPath&&) = default;
    BeamPath &operator=(const BeamPath&) = default;
    BeamPath &operator=(BeamPath&&) = default;

    int getPathMatIndex() const { return pathmat; }
    int getPathPatIndex() const { return pathpat; }

    const HaplotypeHistory &getPathHistoryMat() const { return happaths[pathmat]->getPathHistory(); }
    const HaplotypeHistory &getPathHistoryPat() const { return happaths[pathpat]->getPathHistory(); }
    fp_type getPathProbability() const { return prob; }
    fp_type getPathWeight() const { return weight; }
    int getPathNumErrors() const { return nErrPath; }
    fp_type getPathErrorProb() const { return pErrPath; }
#ifdef DEBUG_TARGET
    int getPathNumForgottenErrors() const { return forgottenerr; }
#endif

    // returns the haplotype pair at the position 'index' sites before the latest position
    pair<Haplotype,Haplotype> getHapsAt(size_t index) const { return make_pair(happaths[pathmat]->getHapAt(index), happaths[pathpat]->getHapAt(index)); }

    // returns the path probability of the selected extension
    fp_type getExtensionProbability(int extensionIndex) const { return extensionProbs[extensionIndex]; }

    // merges the other path to this path not considering extensions,
    // only current probabilities, weight, history etc.
    void merge(BeamPath &other);

    // applies a prepared extension on this BeamPath object and clears the remainder
    // also applies an extension scaling, if provided != 0
    // (also includes the error adjustment (scaling), i.e. if 'deleteErr' is true, one error in the path is removed ('forgotten'))
    void applyExtension(int mat, int pat, int extensionIndex, int weightScaleExp, bool deleteErr);

    // marks an extension as deleted
    void markExtensionDeleted(int extensionIndex);

    // returns true if the selected extension is not deleted, i.e. active
    bool isExtensionActive(int extensionIndex) const { return activeExtensions[extensionIndex]; }

    // returns the next highest active extension index, -1 if non is active
    int getActiveExtensionIndex() { return activeExtensionIndex; }
    // returns the next highest active extension index lower than the parameter, -1 if non is active
    int getNextActiveExtensionIndex(int i);

//    // returns the number of active extensions
//    int getActiveExtensionCount() const { return activeExtensionCount; }

    // mark all required extensions in the corresponding happaths
    void markRequiredExtensions(const Constraint &constr);

    // prepares all required extensions and calculates their probabilities
    // the extensions are:
    // - the four possible haplotype pair extensions (00,01,10,11)
    // the method requires the corresponding haplotype paths already to be extended
    // and takes the target genotype, the current constraint and
    // a global beampath scaling factor (that should be used for all beampaths at this site) as argument.
//    // returns the maximum probability of all extensions (including weight)
    // returns the maximum probability of all extensions (excluding weight)
    fp_type extend(Genotype targetGeno, const Constraint &constr, int32_t bpscaleexp);

    bool operator<(const BeamPath &other) const;

//    // DEBUG
//    HaplotypePathProbability getMatPathProb() { return happaths[pathmat]->getPathProbability(); }
//    HaplotypePathProbability getPatPathProb() { return happaths[pathpat]->getPathProbability(); }
//    HaplotypePathProbability getExtMatPathProb(int extensionIndex) { return happaths[pathmat]->getExtensionProbability(extensionIndex/2); }
//    HaplotypePathProbability getExtPatPathProb(int extensionIndex) { return happaths[pathpat]->getExtensionProbability(extensionIndex%2); }
//    bool isExtMatPathActive(int extensionIndex) { return happaths[pathmat]->isExtensionActive(extensionIndex/2); }
//    bool isExtPatPathActive(int extensionIndex) { return happaths[pathpat]->isExtensionActive(extensionIndex%2); }
//    // __DEBUG

private:

    // after an extension, this method updates the number of errors and the error probability resulting from the extension
    // (also includes the error adjustment (scaling), i.e. if 'deleteErr' is true, one error in the path is removed ('forgotten'))
    inline void updateErrors(int extensionIndex, bool deleteErr);

    inline pair<Haplotype,Haplotype> getHapsFromExtensionIndex(int extensionIndex) const;
    inline bool inconsistent(int extIdentifier, Genotype targetGeno) const;
    inline int errorsToAdd(int extInd, Genotype targetGeno) const;

    inline void extendSingle(int id, fp_type p1, fp_type p2, fp_type scale, bool incon, fp_type &maxP);

    inline void updateActiveExtensionIndex();

    vector<HaplotypePath*> &happaths;
    int pathmat;
    int pathpat;

    fp_type prob;   // probability of the hap pair at the beam front
    fp_type weight; // weight for the hap pair at the beam front

    const fp_type pErr; // general error probability
    int    nErrPath; // current number of errors in path
    fp_type pErrPath; // current error probability factor according to current number of errors in path (pow(pErr,nErr))
    fp_type pErrPathp1; // current error probability factor with one more error
#ifdef DEBUG_TARGET
    int    forgottenerr; // current number of forgotten errors in path (by common error adjustment)
#endif

    array<fp_type, 4> extensionProbs;   // probabilities of all sixteen (4x4) possible extensions (calculated in extend())
    array<bool, 4> activeExtensions; // which extensions are valid, and which deleted
    int activeExtensionIndex; // the highest index of all active extensions (-1 if non active)
//    int activeExtensionCount; // the number of active extensions

    Genotype extTargetGeno; // the target genotype at the extension site

};


#endif /* BEAMPATH_H_ */
