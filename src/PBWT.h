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

#ifndef PBWT_H_
#define PBWT_H_

#include <vector>
#include <iostream>

#include "Datatypes.h"
#include "RingBuffer.h"
#include "PBWTInterval.h"

#include "MyMalloc.h"

using namespace std;

class PBWT {
public:

    // constructor that creates a PBWT from a vector of boolean vectors (natural order)
    // all fields had to be preinitialized with zeros!!!
    PBWT(const RingBuffer<BooleanVector> *dataraw, vector<int> *gCount0, size_t ntarget = -1);
    PBWT(const RingBuffer<BooleanVector> *dataraw, bool refPBWT, size_t ntarget = -1);

    // constructor that takes a preliminary created compressed PBWT as input, that is generated without count0 information (also not in the raw data!)
    PBWT(uint32_t* comprPBWTraw, size_t K, size_t Kwords, size_t M, size_t ntarget = -1);

    // constructor that takes a preliminary created compressed PBWT as input
    PBWT(uint32_t* comprPBWT, vector<int> *gCount0, size_t K, size_t Kwords, size_t M, size_t ntarget = -1);

    ~PBWT() {
        MyMalloc::free(comprPBWT);
        delete gCount0;
    }

    // the forward PBWT will be set to reverse, the current cursor is M and all mappings are expected to be backwards
    void switchToReverse();
    bool isReverse() const { return reversePBWT; }

    void advanceTo(size_t m);

    size_t getK() const { return K; }
    size_t getM() const { return M; }
    fp_type getKInv() const { return Kinv; }

    // returns the number of reference sequences that end with 0 at the provided site
    int getCount0(size_t m) const { return (*gCount0)[m]; }

    // maps an interval from position 'm-1' to all intervals at position 'm',
    // assuming that the sequences indicated by the provided interval are extended with 0 and 1 respectively.
    // 'from' is the provided interval from previous cursor position,
    // 'to0' holds the mapped interval for 0-extension after return,
    // 'to1' holds the mapped interval for 1-extension after return.
    void mapInterval(size_t m, const PBWTInterval &from, PBWTInterval &to0, PBWTInterval &to1) const;

    // returns the PBWT sort order (i.e. the global absolute permutation)...
    // ...before(!) the provided site if this is a forward PBWT or
    // ...at(!) the provided site if this is a backward PBWT
    // NOTE: only possible if storeAbsPerm == true
    const vector<int> &getSortOrder(size_t m) const { return gAbsPerm[m]; }

private:

    // initialization called from constructor
    void init();

    // generate compressed PBWT
    void buildComprPBWT(bool reverse);

    // maps an interval given the part of the compressed PBWT at the desired site
    inline void mapInterval(const uint32_t *comprPBWTatM, int cnt0, const PBWTInterval &from, PBWTInterval &to0, PBWTInterval &to1) const;

    // returns the relative permutation of element "from" as encoded in the compressed PBWT, gcnt0 is the global count0 at this position
    inline int decodeRelPerm(const uint32_t *comprPBWTatM, int gcnt0, int from) const;

    const size_t K; // number of conditioning haps
    const size_t Kwords; // number of 32 bit words used for one site (i.e. at least 2*ceil(K/32), may include padding): used for offsets in the compressed PBWT when jumping to a specific site
    const fp_type Kinv; // inverted number of conditioning haps (1/K)
    const size_t M; // PBWT width (2*number of sites in the condensed reference + 1)
    const RingBuffer<BooleanVector> *dataraw; // condensed reference information of conditioning haps at split sites (transposed)
    uint32_t *comprPBWT;
    bool reversePBWT;
    size_t curr_m; // current cursor

    // relative permutations (global)
    vector<int> *gCount0; // number of reference sequences that end with 0 for all sites

    // absolute permutations (global) - only required (and completetly initialized) for imputation, otherwise the size is 2 for temporary use
    vector<vector<int>> gAbsPerm;
    bool refPBWT;
    bool haveCount0;

    size_t ntarget;

//    // DEBUG
//    void dump(bool offsets);
};

ostream &operator<<(ostream &out, const PBWTInterval &i);

#endif /* PBWT_H_ */
