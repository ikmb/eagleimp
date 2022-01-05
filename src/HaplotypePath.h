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

#ifndef HAPLOTYPEPATH_H_
#define HAPLOTYPEPATH_H_

#include <utility>
#include <memory>
#include <array>
#include <vector>
#include <atomic>

#include "History.h"
#include "Datatypes.h"
#include "HaplotypeHistory.h"
#include "PBWT.h"

#ifdef DEBUG_MAPPINGS_TO_STDERR
struct DbgMapContainer {
    unsigned long long history;
    int history_length;
    int mapping[2];
    bool operator<(const DbgMapContainer &other) const;
    bool operator==(const DbgMapContainer &other) const;
};
#endif

struct HaplotypePathProbability {
    fp_type p;
    int32_t scaleexp;
    HaplotypePathProbability(fp_type p_, int scaleexp_) : p(p_), scaleexp(scaleexp_) {}
};

class HaplotypePath {

    // absolute exponent value must be greater than this before path scaling is applied
    static const int MINSCALEOFFSET;

public:
    // creates a path containing one site with the given history size andpro-forma hap before(!) the first site
    HaplotypePath(size_t historySize, PBWT &pbwt, Haplotype hap);
    // special constructor that copies the given HaplotypePath, but applies the selected (and prepared) extension right away
    // the selected extension in the other object will become invalid afterwards
    HaplotypePath(const HaplotypePath& other, int extensionIndex);
    // default copy + move constructors, copy + move assignments
    HaplotypePath(const HaplotypePath&) = default;
    HaplotypePath(HaplotypePath&&) = default;
    HaplotypePath &operator=(const HaplotypePath&) = default;
    HaplotypePath &operator=(HaplotypePath&&) = default;

    const HaplotypeHistory &getPathHistory() const { return path; }
    // return probability
    HaplotypePathProbability getPathProbability() const { return prob; }

    // returns the haplotype at the position 'index' sites before the latest position
    Haplotype getHapAt(size_t index) const { return path.at(index); }

    // returns the path probability of the selected extension
    HaplotypePathProbability getExtensionProbability(int extensionIndex) const { return HaplotypePathProbability(extensionProbs[extensionIndex], extscaleexp); }
    int getExtScaleExp() const { return extscaleexp; }

    // applies a prepared extension on this HaplotypePath object and clears the remainder
    void applyExtension(int extensionIndex);

    // marks an extension as deleted
    void markExtensionDeleted(int extensionIndex);

    // returns true if the selected extension is not deleted, i.e. active
    bool isExtensionActive(int extensionIndex) const { return activeExtensions[extensionIndex]; }

    bool allExtensionsActive() const { return activeExtensions[0] && activeExtensions[1]; }
    bool allExtensionsInactive() const { return !(activeExtensions[0] || activeExtensions[1]); }

    // returns true if it is required to extend this happath
    bool isToBeExtended() const { return toextend[0] || toextend[1]; }

    // marks the extensions that have to be calculated later using the provided constraint
    void markRequiredExtensions(const Constraint &constr);

    void incExtReferenceCnt(int extensionIndex) { ext_reference_cnt[extensionIndex]++; }

    // prepares all required extensions (previously marked) and calculates their probabilities
    // uses the provided vector of recombination probabilities (latest is the actual one for the site 'cursor' (refers to extension site, split sites domain))
    // if copyPrevHistSites is >0, that many sites of the history of prevPath are copied instead of newly calculated.
    // returns the scaling exponent used for the extensions.
    int32_t extend(const vector<fp_type> &prec, size_t cursor,
            const HaplotypePath *prevPath
#if defined DEBUG_TARGET && defined DEBUG_PRINT_PBWT_MAPPINGS
            , const Constraint &constr
#endif
#ifdef DEBUG_MAPPINGS_TO_STDERR
            , vector<DbgMapContainer> &dbgmappings
#endif
            );

    // returns true if all stored interval mappings are equal
    bool isMergeable(const HaplotypePath &other) const;

//    bool operator<(const HaplotypePath &other) const;

#ifdef DEBUG_PROFILE_PBWT_MAPPINGS
    static array<atomic<size_t>,32> dbg_mapstop;
    static array<atomic<size_t>,32> dbg_mapone;
    static atomic<size_t> dbg_mapstop_gt;
    static atomic<size_t> dbg_mapone_gt;
#endif

private:
//    using pbwtInts_type = vector<Interval>;
//    using pbwtIntsPtr_type = shared_ptr<pbwtInts_type>;

    // adjusts scaling parameters according to current path probability
    void updateScaling();

    // calculate P(h_{1:m-1}0) and P(h_{1:m-1}1)
    // the results are returned in p0 for zero extension, and p1 for one extension respectively
    void calculateHaplotypeProbabilities(
            const vector<fp_type> &pRec,         // recombination probabilities
            size_t cursor,                      // current site 'm'
            bool calcP0, // if desired, one could omit the calculation of either probability to save runtime
            bool calcP1, // by setting the corresponding flag to 'false'
            fp_type &p0,    // the result for zero extension via a consistent segment
            fp_type &p1,    // the result for one extension via a consistent segment
            vector<PBWTInterval> &extPbwtIntPtr0,    // reference to vector that will hold the new intervals afterwards for 0 extension (preliminary contents will be overwritten!)
            vector<PBWTInterval> &extPbwtIntPtr1,    // reference to vector that will hold the new intervals afterwards for 1 extension (preliminary contents will be overwritten!)
            size_t copyPrevHistSites, // this many sites are copied from intervals below (to save the new calculation)
            const HaplotypePath *prevPath
#ifdef DEBUG_MAPPINGS_TO_STDERR
            , vector<DbgMapContainer> &dbgmappings
#endif
            ) const;

    inline Haplotype getHapFromExtensionIndex(int extensionIndex) const { return extensionIndex == 0 ? Haplotype::Ref : Haplotype::Alt; }

    size_t historySize;
    // the haplotype history of the path
    HaplotypeHistory path;
    HaplotypePathProbability prob;   // probability of the hap path
    // probabilities of hap path at each position in history
    History<fp_type> probHist;

    PBWT &pbwt; // PBWT structure used for frequency lookups
    // intervals where subsequences of the path are found in the PBWT structure
    vector<PBWTInterval> pbwtInt;

    // probabilities of all extensions calculated in extend()
    array<fp_type, 2> extensionProbs;
    int extscaleexp; // scaling exponent of the extensions
    // PBWT intervals of the extensions
    array<vector<PBWTInterval>, 2> extensionPbwtInt;

    bool activeExtensions[2]; // which extensions are valid, and which deleted/inactive.
    bool toextend[2];         // which extensions should be processed
    int ext_reference_cnt[2]; // keeps track of how many beam path extensions reference to this haplotype path 0/1-extension

#ifdef DEBUG_MAPPINGS_TO_STDERR
    uint64_t dbgmap_currhistory = 0ull;
    uint64_t dbg_path2uint64(const HaplotypeHistory&) const;
#endif
};

//using HaplotypePathPtr = shared_ptr<HaplotypePath>;

#endif /* HAPLOTYPEPATH_H_ */
