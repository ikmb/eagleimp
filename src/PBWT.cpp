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
#include <vector>

#include "utils.h"

#include "PBWT.h"

#ifdef STOPWATCH
#include "Stopwatch.h"
#endif

PBWT::PBWT(const RingBuffer<BooleanVector> *dataraw_, vector<int> *gCount0_, size_t ntarget_)
    : K((*dataraw_)[0].size()),
//      Kwords(4*roundToMultiple((size_t)K, UNITWORDS * sizeof(BooleanVector::data_type) * 8) / (8 * sizeof(BooleanVector::data_type))),
      Kwords(2*divideRounded(K,(size_t)32)),
      Kinv(1.0/(fp_type)K),
      M(dataraw_->size()),
      dataraw(dataraw_),
      reversePBWT(false),
      curr_m(-1),
      refPBWT(false),
      haveCount0(true),
      ntarget(ntarget_)
{
    gCount0 = gCount0_;
    init();

//    // DEBUG
//    dump(ntarget, true);
}

PBWT::PBWT(const RingBuffer<BooleanVector> *dataraw_, bool refPBWT_, size_t ntarget_)
    : K((*dataraw_)[0].size()),
      Kwords(2*divideRounded(K,(size_t)32)),
      Kinv(1.0/(fp_type)K),
      M(dataraw_->size()),
      dataraw(dataraw_),
      reversePBWT(false),
      curr_m(-1),
      refPBWT(refPBWT_),
      haveCount0(false),
      ntarget(ntarget_)
{
    gCount0 = new vector<int>(M,0);
    init();

//    // DEBUG
//    dump(ntarget, true);
}

PBWT::PBWT(uint32_t* comprPBWTraw, size_t K_, size_t Kwords_, size_t M_, size_t ntarget_)
    : K(K_),
      Kwords(2*Kwords_), // will be the size after filling with cnt0 offsets
      Kinv(1.0/(fp_type)K),
      M(M_),
      dataraw(NULL),
      reversePBWT(false),
      curr_m(-1),
      refPBWT(false),
      haveCount0(true),
      ntarget(ntarget_)
{
    gCount0 = new vector<int>(M,0);
    comprPBWT = (uint32_t*) MyMalloc::malloc(2*Kwords*M*(32/8), string("comprPBWT_t")+to_string(ntarget_)); // generate space also for the count0 offsets and fwd/bck PBWT

    uint32_t *src = comprPBWTraw;
    uint32_t *dst = comprPBWT;
    for (size_t m = 0; m < 2*M; m++) {
        uint32_t currcnt0 = 0;
        for (size_t k = 0; k < Kwords_; k++) { // Note: iterating over the original number of Kwords!
            currcnt0 += 32-__builtin_popcount(*src);
            *dst++ = *src++; // copy data
            *dst++ = currcnt0; // insert cnt0 offset
        }
    }

    uint32_t *currcnt0ptr = comprPBWT + (2*divideRounded(K,(size_t)32)-1); // points to cnt0 of first site
//            // DEBUG
//            cerr << "cnt0s " << tgt->getIdx() << ": " << tgt->getNSplits() << " " << M << endl; size_t idx = 0;
    for (auto &gcnt0 : *gCount0) {
        gcnt0 = *currcnt0ptr;
        currcnt0ptr += Kwords;
//                // DEBUG
//                cerr << (idx++) << ": " << gcnt0 << endl;
    }

//    // DEBUG
//    dump(ntarget, true);
}

PBWT::PBWT(uint32_t* comprPBWT_, vector<int> *gCount0_, size_t K_, size_t Kwords_, size_t M_, size_t ntarget_)
    : K(K_),
      Kwords(Kwords_),
      Kinv(1.0/(fp_type)K),
      M(M_),
      dataraw(NULL),
      reversePBWT(false),
      curr_m(-1),
      refPBWT(false),
      haveCount0(true),
      ntarget(ntarget_)
{
    gCount0 = gCount0_;
    comprPBWT = comprPBWT_;

//    // DEBUG
//    dump(ntarget, true);
}

void PBWT::init() {
#ifdef STOPWATCH
    Stopwatch swpbwt("PBWT");
#endif

    comprPBWT = (uint32_t*) MyMalloc::malloc((refPBWT ? 1 : 2)*M*Kwords*sizeof(uint32_t), string("comprPBWT_t")+to_string(ntarget)); // space for fwd and bck PBWT (latter only if not for reference)

    // depending on if the absolute permutations are to be stored, reserve memory either for temporary or global usage
    if (refPBWT)
        gAbsPerm = move(vector<vector<int>>(M+1, vector<int>(K,0)));
    else
        gAbsPerm = move(vector<vector<int>>(2, vector<int>(K,0))); // need only "actual" and "new"

    // build forward and backward PBWT
    buildComprPBWT(false); // forward
    haveCount0 = true;
    if (!refPBWT) // TODO if we don't do reverse phasing, we don't need this!
        buildComprPBWT(true); // backward

#ifdef STOPWATCH
    swpbwt.stop();
#endif
}

//// DEBUG dump PBWT
//void PBWT::dump(bool offsets) {
//
//    cout << "DEBUG: dump " << ntarget << endl;
//    ofstream ofs(string("pbwt")+to_string(ntarget));
//    const size_t Krnd = 2*divideRounded(K,32);
//    const size_t Kpad = Kwords - Krnd;
//    ofs << "K: " << K << " Kpbwt: " << Kwords << " Kpad: " << Kpad << endl;
//    ofs << "Mspl: " << M/2 << " Mpbwt: " << M << endl;
//    uint32_t* dptr = comprPBWT;
//    for (int fwdbck = 0; fwdbck < 2; fwdbck++) {
//        if (fwdbck == 0)
//            ofs << "fwd:" << endl;
//        else
//            ofs << "bck:" << endl;
//        for (size_t midx = 0; midx < M; midx++) {
//            size_t m = fwdbck ? M-midx-1 : midx;
//            ofs << dec << m << ":";
//            int count0 = (*gCount0)[m];
//            int count0x = dptr[Krnd-1];
//            ofs << " count0: " << count0 << (count0 == count0x ? " ok" : " XXX");
////            for (size_t k = 0; k < Kwords/2; k++) {
//            for (size_t k = 0; k < Krnd/2; k++) {
//                if (k%8 == 0)
//                    ofs << endl;
//                ofs << " " << hex << setw(8) << setfill('0') << *dptr++; // data
//                if (offsets)
//                    ofs << " " << dec << setw(5) << *dptr++; // offset
//                else
//                    ofs << " " << hex << setw(8) << setfill('0') << *dptr++; // print offset field as data
//            }
//            for (size_t k = Krnd; k < Kwords; k++, dptr++);
//            ofs << endl;
//        }
//    }
//    ofs.close();
//
//};
//// __DEBUG


void PBWT::switchToReverse() {
    // switch to reverse
    reversePBWT = true;
    curr_m = M;
}

void PBWT::advanceTo(size_t m) {
    curr_m = m;
}

// NOTE: takes advantage of pre-calculated count0 if possible
void PBWT::buildComprPBWT(bool reverse) {
    uint32_t *compr_ptr = comprPBWT;
    if (reverse)
        compr_ptr += M*Kwords;

    // the first absPerm is the identity
    size_t idx = reverse && refPBWT ? M : 0;
    for (size_t i = 0; i < K; i++)
        gAbsPerm[idx][i] = i;

    for (size_t midx = 0; midx < M; midx++) {

        uint32_t *curr_compr_ptr = compr_ptr;

        size_t m = reverse ? M-1-midx : midx;
        const BooleanVector &curr_data = (*dataraw)[m];

        // calculate count0:
        // we can determine the number of zeroes in advance, which is advantageous.
        int cnt0 = haveCount0 ? (*gCount0)[m] : curr_data.getFalseCount_withPreInit();

        size_t abspermidx = refPBWT ? (reverse ? m+1 : m) : 0;
        size_t nabspermidx = refPBWT ? (reverse ? m : m+1) : 1;
        size_t p0next = 0;
        size_t p1next = cnt0;
        auto it = gAbsPerm[abspermidx].cbegin();
        auto &absPermNew = gAbsPerm[nabspermidx];
        uint32_t curr_permdata = 0; // current vector in permuted order
        uint32_t curr_mask = 1;
        for (size_t i = 0; i < K; i++, ++it) {
            if (!curr_data[*it]) { // curr sample is 0

                // insert index of this 0-sequence in 0 part of vector
                absPermNew[p0next] = *it;
                p0next++;

            } else { // curr sample is 1

                curr_permdata |= curr_mask; // store 1 in current perm vector

                // insert index of this 1-sequence in 1 part of vector
                absPermNew[p1next] = *it;
                p1next++;
            }
            curr_mask <<= 1; // select next bit
            if (curr_mask == 0) { // the '1' has been shifted out of the vector
                *curr_compr_ptr = curr_permdata; // store permuted vector
                curr_compr_ptr++; // now points to offset field
                *curr_compr_ptr = p0next; // store current cnt0 (which is p0next) as offset
                curr_compr_ptr++; // now points to next data field
                curr_permdata = 0; // reset permutation data
                curr_mask = 1; // reset mask
            }
        }
        if (K%32) { // need to write the last incomplete word
            curr_mask = ~(curr_mask-1); // set all remaining bits to one in order to insert 1-padding
            *curr_compr_ptr = curr_permdata | curr_mask; // insert 1-padding (necessary for interval mapping!)
            curr_compr_ptr++;
            *curr_compr_ptr = p0next; // same as cnt0 now
            curr_compr_ptr++;
        }

        // if we need to add some padding to the compressed PBWT, we must do this here
        compr_ptr += Kwords; // independent if we need padding or not, or which part we're on, this'll always do!

        // store count0 (if we still need it)
        if (!haveCount0)
            (*gCount0)[m] = p0next;

        // Only if absPerm is not stored globally:
        // The just calculated absPermNew will be the absPerm of the next step. So, we move both vectors that we don't need to re-allocate memory.
        if (!refPBWT) {
            vector<int> tmp(move(gAbsPerm[0]));
            gAbsPerm[0] = move(gAbsPerm[1]);
            gAbsPerm[1] = move(tmp);
        }
    }

}


void PBWT::mapInterval(size_t m, const PBWTInterval &from, PBWTInterval &to0, PBWTInterval &to1) const {

    mapInterval(&comprPBWT[reversePBWT ? (2*M-1-m)*Kwords : m*Kwords], (*gCount0)[m], from, to0, to1);

#if defined DEBUG_TARGET and defined DEBUG_PRINT_PBWT_MAPPINGS
    cout << "map: " << from << " -> "
            << to0 << " + " << to1 << " cnt0: " << gCount0[m] << endl;
#endif

}


inline void PBWT::mapInterval(const uint32_t *comprPBWTatM, int cnt0, const PBWTInterval &from, PBWTInterval &to0, PBWTInterval &to1) const {
    // lookup from-sites and calculate to-sites

    int fl = from.getLeft();
    int fr = from.getRight();

    // default
    to0.setEmpty();
    to1.setEmpty();

    if (from.isNotEmpty()) {
        if (from.isSingle()) { // if the interval contains only one element, we can speed-up the process a bit

            int i = decodeRelPerm(comprPBWTatM, cnt0, fl); // where does the element map to?
            if (i < cnt0) { // corresponding sequence has '0' at cursor
                // so this is a 0-extension
                to0.set(i, i);
                // the 1-extension will be empty
            } else { // corresponding sequence has '1' at cursor
                // so this is a 1-extension
                to1.set(i, i);
                // the 0-extension will be empty
            }

        } else { // interval contains more than one element

#if defined DEBUG_MAPPINGS_TO_STDERR || (defined DEBUG_TARGET && defined  DEBUG_PRINT_PBWT_MAPPINGS)
            bool swappedl = false, swappedr = false;
#endif

            // map into haplotype-permutation

            // left border
            int l0 = decodeRelPerm(comprPBWTatM, cnt0, fl); // mapping of left border
            int l1 = fl - l0 + cnt0; // alternative mapping if bit at cursor was inverted
            if (l0 >= cnt0) { // sequence of left border has '1' at current cursor
                swap(l0, l1); // direct mapping is left border of '1'-extension, so swap borders
#if defined DEBUG_MAPPINGS_TO_STDERR || (defined DEBUG_TARGET && defined  DEBUG_PRINT_PBWT_MAPPINGS)
                swappedl = true;
#endif
            }
            // else sequence of left border has '0' at current cursor, and everything's fine

            // right border
            int r0 = decodeRelPerm(comprPBWTatM, cnt0, fr); // mapping of right border
            int r1 = fr - r0 + cnt0 - 1; // alternative mapping if bit at cursor was inverted
            if (r0 >= cnt0) { // sequence of right border has '1' at current cursor
                swap(r0, r1);  // direct mapping is right border of '1'-extension, so swap borders
#if defined DEBUG_MAPPINGS_TO_STDERR || (defined DEBUG_TARGET && defined  DEBUG_PRINT_PBWT_MAPPINGS)
                swappedr = true;
#endif
            }
            // else sequence of right border has '0' at current cursor, and everything's fine

            // return values
            to0.set(l0, r0);
            to1.set(l1, r1);
#if defined DEBUG_MAPPINGS_TO_STDERR || (defined DEBUG_TARGET && defined  DEBUG_PRINT_PBWT_MAPPINGS)
            to0.altl = swappedl;
            to0.altr = swappedr;
            to1.altl = not swappedl;
            to1.altr = not swappedr;
#endif
        } // END more than one element
    } // END is not empty
}

inline int PBWT::decodeRelPerm(const uint32_t *comprPBWTatM, int gcnt0, int from) const {
    uint32_t off = comprPBWTatM[2*(from/32)+1]; // offset: count0 including the current field
    uint32_t field = comprPBWTatM[2*(from/32)]; // current data field
    // to restore cnt0 until current position, need to subtract 0's following the current position including current position
    uint32_t mask = (1 << (from%32)); // marks the actual bit
    uint32_t bit = field & mask; // extract current value
    mask--; // all bits up to current position (excluding curr pos) are 1
    field |= mask; // current field with all lower bits set to 1 -> all 0's in field need to be subtracted from offset to get current cnt0
    uint32_t fcnt0 = off - __builtin_popcount(~field); // subtract 0's by counting 1's in inverse field
    if (bit) { // current value is 1
        // need to calculate cnt1 up to here: must be all haps minus cnt0
        int fcnt1 = from - (int)fcnt0;
        return gcnt0 + fcnt1;
    } else // current value is 0
        return (int)fcnt0; // cnt0 up to here is the new relative permutation index
}

ostream &operator<<(ostream &out, const PBWTInterval &i) {
    if (i.isEmpty())
        out << "[empty]";
    else {
#if defined DEBUG_MAPPINGS_TO_STDERR || (defined DEBUG_TARGET && defined  DEBUG_PRINT_PBWT_MAPPINGS)
        out << "[" << i.getLeft();
        if (i.altl)
            out << "*";
        out << ", " << i.getRight();
        if (i.altr)
            out << "*";
        out << "](" << i.size() << ")";
#else
        out << "[" << i.getLeft() << ", " << i.getRight() << "](" << i.size() << ")";
#endif
    }
    return out;
}

