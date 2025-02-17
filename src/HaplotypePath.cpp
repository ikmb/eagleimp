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

#include "HaplotypePath.h"

/*static*/ const int HaplotypePath::MINSCALEOFFSET = 20;
#ifdef DEBUG_PROFILE_PBWT_MAPPINGS
    /* static */ array<atomic<size_t>,32> HaplotypePath::dbg_mapstop;
    /* static */ array<atomic<size_t>,32> HaplotypePath::dbg_mapone;
    /* static */ atomic<size_t> HaplotypePath::dbg_mapstop_gt = 0;
    /* static */ atomic<size_t> HaplotypePath::dbg_mapone_gt = 0;
#endif

HaplotypePath::HaplotypePath(size_t historySize_, PBWT &pbwt_, Haplotype hap)
 : historySize(historySize_),
   path(historySize_),
   prob(1.0, 0), // empty path has probability 1.0
   probHist(historySize_),
   pbwt(pbwt_),
   activeExtensions{false,false},
   toextend{false,false},
   ext_reference_cnt{false,false}
{
    // haplotype history containing the given haplotype as pro-forma path before(!) the first site
    path.push_new(hap);

    // empty path has 1.0 probability
    probHist.push_new(1.0);

    // clear extension fields for new four possible extensions
    extensionProbs.fill(0.0);
    extscaleexp = 0;
}

HaplotypePath::HaplotypePath(const HaplotypePath &other, int extensionIndex)
 : historySize(other.historySize),
   path(other.path), //errors(other.errors),
   prob(other.extensionProbs[extensionIndex], other.extscaleexp), // set the probability to that of the selected extension
   probHist(other.probHist),
   pbwt(other.pbwt),
   pbwtInt(move(other.extensionPbwtInt[extensionIndex])), // move! extension is invalid afterwards!
   activeExtensions{false,false},
   toextend{false,false},
   ext_reference_cnt{false,false}
{
    // set the haps of the selected extension to the end of the path
    path.push_new(getHapFromExtensionIndex(extensionIndex));

    // push the selected extension probabilities to the end of the history
    probHist.push_new(other.extensionProbs[extensionIndex]);

    // clear extension fields for new possible extensions
    extensionProbs.fill(0.0);
    extscaleexp = 0;
}

void HaplotypePath::applyExtension(int extensionIndex) {

    // set the haps of the selected extension to the end of the path
    path.push_new(getHapFromExtensionIndex(extensionIndex));

    // set the probability to that of the selected extension
    prob = getExtensionProbability(extensionIndex);

    // push the selected extension probabilities to the end of the history
    probHist.push_new(extensionProbs[extensionIndex]);

    // take over the intervals from the corresponding extension (move! extension is invalid afterwards!)
    pbwtInt = move(extensionPbwtInt[extensionIndex]);

    // clear extension fields for new possible extensions
    extensionProbs.fill(0.0);
    extscaleexp = 0;

    activeExtensions[0] = false;
    activeExtensions[1] = false;
    toextend[0] = false;
    toextend[1] = false;
    ext_reference_cnt[0] = 0;
    ext_reference_cnt[1] = 0;
}

void HaplotypePath::markExtensionDeleted(int extensionIndex) {
    ext_reference_cnt[extensionIndex]--;
    if (ext_reference_cnt[extensionIndex] == 0)
        activeExtensions[extensionIndex] = false;
}

#ifdef DEBUG_MAPPINGS_TO_STDERR
bool DbgMapContainer::operator<(const DbgMapContainer &other) const {
    if (history_length < other.history_length)
        return true;
    else if(history_length == other.history_length && history < other.history)
        return true;
    else
        return false;
}
bool DbgMapContainer::operator==(const DbgMapContainer &other) const {
    return history_length == other.history_length && history == other.history; // the mappings should be equal as well
}
uint64_t HaplotypePath::dbg_path2uint64(const HaplotypeHistory& dbgpath) const {
    uint64_t ret = 0ull;
    for (size_t i = 0; i < min(dbgpath.size(), (size_t)64); i++) {
        Haplotype h = dbgpath.at(i);
        if (h == Haplotype::Alt) {
            ret |= 1ull << i;
        }
    }
    return ret;
}
#endif

void HaplotypePath::markRequiredExtensions(const Constraint &constr) {

    // DO NOT MARK ANYTHING WITH FALSE HERE SINCE IT MAY ALREADY BE MARKED BY ANOTHER PATH!!!
    if (constr.type == Constraint::Type::NoConstr) {
        toextend[0] = true;
        toextend[1] = true;
    } else {
        if ((constr.type == Constraint::Type::HetKeep && path[constr.refidx] == Haplotype::Ref)
                || (constr.type == Constraint::Type::HetSwitch && path[constr.refidx] == Haplotype::Alt)
                || (constr.type == Constraint::Type::HomRef)) { // 0
            toextend[0] = true;
        } else if ((constr.type == Constraint::Type::HetKeep && path[constr.refidx] == Haplotype::Alt)
                || (constr.type == Constraint::Type::HetSwitch && path[constr.refidx] == Haplotype::Ref)
                || (constr.type == Constraint::Type::HomAlt)) { // 1
            toextend[1] = true;
        }
    }
}

// returns the new scaling exponent for the extensions
int32_t HaplotypePath::extend(const vector<fp_type> &prec, size_t cursor,
        const HaplotypePath *prevPath
#if defined DEBUG_TARGET && defined DEBUG_PRINT_PBWT_MAPPINGS
        , const Constraint &constr
#endif
#ifdef DEBUG_MAPPINGS_TO_STDERR
        , vector<DbgMapContainer> &dbgmappings
#endif
        ) {
    // adjust scaling before calculating the new probabilities
    updateScaling();

    // calculate all probabilities for the two possible extensions

    // determine number of equal haps from previous path
    size_t copyPrevHistSites = 0ull;
    if (prevPath) { // not null
        size_t numequal = path.numEqual(prevPath->path);
        copyPrevHistSites = min(historySize, 1 + numequal);
    }

#if defined DEBUG_TARGET && defined DEBUG_PRINT_PBWT_MAPPINGS
    cout << "  Constraint: " << constr << " copy: " << copyPrevHistSites << endl;
    cout << "  path: ";
    path.printPath(cout, -1);
    cout << endl;
#endif

    // new PBWT intervals
    vector<PBWTInterval> extint0;
    vector<PBWTInterval> extint1;

    extint0.reserve(historySize);
    extint1.reserve(historySize);

    // select path probabilities that have to be calculated according to previous markings
    bool calcp0 = toextend[0], calcp1 = toextend[1];

    // calculate path probabilities if extended with 0/1
#ifdef DEBUG_MAPPINGS_TO_STDERR
    dbgmap_currhistory = dbg_path2uint64(path);
#endif
    fp_type p0, p1;
    calculateHaplotypeProbabilities(prec, cursor, calcp0, calcp1,
            p0,
            p1,
            extint0,
            extint1,
            copyPrevHistSites,
            prevPath
#ifdef DEBUG_MAPPINGS_TO_STDERR
            , dbgmappings
#endif
            );

    extensionProbs[0] = p0;
    extensionProbs[1] = p1;

    extensionPbwtInt[0] = move(extint0);
    extensionPbwtInt[1] = move(extint1);

    activeExtensions[0] = toextend[0];
    activeExtensions[1] = toextend[1];
//    toextend[0] = false; not required
//    toextend[1] = false;

    return extscaleexp;
}

void HaplotypePath::updateScaling() {
    int exp = ilogb(prob.p);
    int newscale = prob.scaleexp; // default is the old scaling value
#ifdef DEBUG_TARGET
    cout << "  usc: " << newscale << " shift: " << exp;
#endif
    if (abs(exp) >= MINSCALEOFFSET) { // we need a new scaling!
        // scaling is adjusted by the exponent
        newscale += exp;
#ifdef DEBUG_TARGET
        cout << " updated: " << newscale << endl;
#endif
        // scaling is achieved by applying the new adjustment to all probs in the history
        // (to save runtime, this is done only to those probs that will be used in further calculations,
        // i.e. the ones where a mapped interval exists for)
        // if the exponent is positive, we have to divide
//        if (exp > 0) {
//            for (size_t i = 0; i <= pbwtInt->size(); i++) { // need to iterate over all intervals + 1 !!
//                probHist[i] /= (fp_type)(1ull << exp);
//            }
//        } else { // negative -> multiply
        // exponent is always negative because scaling is removed from the prec vector
//        assert(exp<0);
            for (size_t i = 0; i <= pbwtInt.size(); i++) { // need to iterate over all intervals + 1 !!
//                probHist[i] *= (fp_type)(1ull << -exp);
                probHist[i] *= pow(2.0,-exp);
            }
//        }
    }
#ifdef DEBUG_TARGET
    else
        cout << " not updated" << endl;
#endif
    extscaleexp = newscale;
}

void HaplotypePath::calculateHaplotypeProbabilities(
        const vector<fp_type> &pRec,
        size_t cursor,
        bool calcP0, bool calcP1,
        fp_type &p0_,
        fp_type &p1_,
        vector<PBWTInterval> &extPbwtInts0,
        vector<PBWTInterval> &extPbwtInts1,
        size_t copyPrevHistSites,
        const HaplotypePath* prevPath
#ifdef DEBUG_MAPPINGS_TO_STDERR
        , vector<DbgMapContainer>& dbgmappings
#endif
        ) const {


    // map all intervals and calculate the probabilities:

    fp_type p0 = 0.0;
    fp_type p1 = 0.0;
    vector<fp_type> tmp(historySize);

    // precalculate all tmp values
    for (size_t i = 0; i < min(historySize, pbwtInt.size()+1); i++)
        tmp[i] = probHist[i] * pRec[i];

    // start with the most recent site (special, since we do not need to map a previous PBWT interval)
    {

        // to determine f(h_n:n) do only mapping for 0 and 1.
        // for f(h_(n-x):n) with x>=1 map only consistent parts (which can be zero -> no mapping then)

        // the first intervals take all 0's and 1's, no inconsistencies to check for
        // path begin is the same as path end
        int cnt0 = pbwt.getCount0(2*cursor+1); // the odd sites mark the reference sites, evens are incon segments
        if (cnt0 > 0) {
            if (calcP0)
                p0 = cnt0 * tmp[0];
            extPbwtInts0.emplace_back(0, cnt0-1);
        }
        if ((size_t)cnt0 < pbwt.getK()) {
            if (calcP1)
                p1 = (pbwt.getK()-cnt0) * tmp[0];
            extPbwtInts1.emplace_back(cnt0, pbwt.getK()-1);
        }

#ifdef DEBUG_PROFILE_PBWT_MAPPINGS
        if (cnt0 == 1)
            dbg_mapone[0]++;
        if (cnt0 == K-1)
            dbg_mapone[0]++;
#endif


#if defined DEBUG_TARGET && defined DEBUG_PRINT_PBWT_MAPPINGS
        cout << "  h: 0 start -> "
                << PBWTInterval(0, cnt0-1) << " + "
                << PBWTInterval(cnt0, pbwt.getK()-1)
                << " p0: " << p0
                << " p1: " << p1
                << endl;
#endif
#ifdef DEBUG_MAPPINGS_TO_STDERR
        DbgMapContainer dbgmc;
        dbgmc.history_length = 0;
        dbgmc.history = 0ull; // no history
        dbgmc.mapping[0] = cnt0;
        dbgmc.mapping[1] = pbwt.getK()-cnt0;
        dbgmappings.push_back(dbgmc);
#endif
    } // end first site

    // flag for preliminary stop
    bool bothempty = false;

    // the history sites that need to be copied
    size_t i = 1;
    for (; !bothempty && i < min(copyPrevHistSites, pbwtInt.size()+1); i++) { // from newest to oldest, but leave out most recent

        bothempty = true;
        // prevPath is not null if we are here (since copyPrevHistSites is non-zero)
        if (prevPath->extensionPbwtInt[0].size() > i) {
            const PBWTInterval &to0 = (prevPath->extensionPbwtInt[0])[i];
            if (calcP0)
                p0 += to0.size() * tmp[i];
            extPbwtInts0.push_back(to0);
            bothempty = false;
        }
        if (prevPath->extensionPbwtInt[1].size() > i) {
            const PBWTInterval &to1 = (prevPath->extensionPbwtInt[1])[i];
            if (calcP1)
                p1 += to1.size() * tmp[i];
            extPbwtInts1.push_back(to1);
            bothempty = false;
        }

#ifdef DEBUG_PROFILE_PBWT_MAPPINGS
        const PBWTInterval &from = pbwtInt[i-1];
        if (from.size() > 1) {
            if (extPbwtInts0.size() == i+1 && extPbwtInts0.back().size() == 1) {
                if (i < dbg_mapone.size())
                    dbg_mapone[i]++;
                else
                    dbg_mapone_gt++;
            }
            if (extPbwtIntPtr1->size() == i+1 && extPbwtIntPtr1->back().size() == 1) {
                if (i < dbg_mapone.size())
                    dbg_mapone[i]++;
                else
                    dbg_mapone_gt++;
            }
        }
#endif

#if defined DEBUG_TARGET && defined DEBUG_PRINT_PBWT_MAPPINGS
        const PBWTInterval &from = pbwtInt[i-1];
        cout << "  h: " << i << " " << from << " -> "
                << (extPbwtInts0.size() == i+1 ? extPbwtInts0.back() : PBWTInterval()) << " + "
                << (extPbwtInts1.size() == i+1 ? extPbwtInts1.back() : PBWTInterval())
                << " p0: " << p0
                << " p1: " << p1
                << endl;
#endif
#ifdef DEBUG_MAPPINGS_TO_STDERR
        // copied intervals are not recorded since they should already be in the list
#endif
    } // end for remaining sites


    // the remaining available history sites that need to be calculated
    for (; !bothempty && i < min(historySize, pbwtInt.size()+1); i++) { // from newest to oldest, but leave out most recent

        // calculate intervals
        const PBWTInterval &from = pbwtInt[i-1]; // need the interval from the same site from the previous iteration

        PBWTInterval to0, to1, toCon, toIncon; // empty

        // check, if previous mapping was equal
        if (i > 1 && pbwtInt[i-2] == from) { // equal mapping, just copy destination
            if (extPbwtInts0.size() == i) // interval exists
                to0.set(extPbwtInts0.back());
            if (extPbwtInts1.size() == i) // interval exists
                to1.set(extPbwtInts1.back());
        } else {
            // PBWT mapping over inconsistency segment (need to know if we are moving forwards or backwards)
            // In the PBWT odd sites mark the reference sites, evens are incon segments
            if (pbwt.isReverse())
                pbwt.mapInterval(2*(cursor+1), from, toCon, toIncon);
            else
                pbwt.mapInterval(2*cursor, from, toCon, toIncon);
            pbwt.mapInterval(2*cursor+1, toCon, to0, to1); // map extensions
        }

#ifdef DEBUG_PROFILE_PBWT_MAPPINGS
        if (from.size() > 1) {
            if (to0.size() == 1) {
                if (i < dbg_mapone.size())
                    dbg_mapone[i]++;
                else
                    dbg_mapone_gt++;
            }
            if (to1.size() == 1) {
                if (i < dbg_mapone.size())
                    dbg_mapone[i]++;
                else
                    dbg_mapone_gt++;
            }
        }
#endif

        bothempty = true;

        // store intervals + calculate probabilities
        if (to0.isNotEmpty()) {
            if (calcP0)
                p0 += to0.size() * tmp[i];
            extPbwtInts0.push_back(to0);
            bothempty = false;
        }
        if (to1.isNotEmpty()) {
            if (calcP1)
                p1 += to1.size() * tmp[i];
            extPbwtInts1.push_back(to1);
            bothempty = false;
        }

#if defined DEBUG_TARGET && defined DEBUG_PRINT_PBWT_MAPPINGS
        cout << "  h: " << i << " ";
//        if (lastfrom == from) {
//           cout << "COPY ";
//        }
           cout << from << " -> "
                << to0 << " + "
                << to1
                << " p0: " << p0
                << " p1: " << p1
                << endl;
#endif
#ifdef DEBUG_MAPPINGS_TO_STDERR
        DbgMapContainer dbgmc;
        dbgmc.history_length = i;
        uint64_t dbgmask = i >= 64 ? ~0ull : ((1ull << i) - 1ull);
        dbgmc.history = dbgmap_currhistory & dbgmask;
        dbgmc.mapping[0] = to0.size();
        dbgmc.mapping[1] = to1.size();
        dbgmappings.push_back(dbgmc);
#endif
    } // end for calculating remaining sites

#ifdef DEBUG_PROFILE_PBWT_MAPPINGS
    if (i-1 < dbg_mapstop.size())
        dbg_mapstop[i-1]++;
    else
        dbg_mapstop_gt++;
#endif

    // returns
    p0_ = p0;
    p1_ = p1;
}

bool HaplotypePath::isMergeable(const HaplotypePath &other) const {
    return path.isEqual(other.path);
}

//bool HaplotypePath::operator<(const HaplotypePath &other) const {
//    return prob > other.prob;
//}
