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

#ifndef PBWTINTERVAL_H_
#define PBWTINTERVAL_H_

#include "Datatypes.h"

// borders inclusive!
class PBWTInterval {
public:
    // empty interval
//    PBWTInterval() : i(0), j(-1), hfrom(Haplotype::Ref), hto(Haplotype::Ref) {}
//    PBWTInterval() : i(0), j(-1), hto(Haplotype::Ref) {}
    PBWTInterval() : i(0), j(-1) {}

    // interval given the borders and haplotype
//    PBWTInterval(int i_, int j_, Haplotype hfrom_, Haplotype hto_) : i(i_), j(j_), hfrom(hfrom_), hto(hto_) {}
//    PBWTInterval(int i_, int j_, Haplotype hto_) : i(i_), j(j_), hto(hto_) {}
    PBWTInterval(int i_, int j_) : i(i_), j(j_) {}

    // get sizes:
    int size()      const { return j - i + 1; }

    // flags
    bool isEmpty()    const { return i > j; }
    bool isNotEmpty() const { return i <= j; }
    bool isSingle()   const { return i == j; }

    void set(int left, int right) {
        i = left;
        j = right;
//        hfrom = hfrom_;
//        hto = hto_;
#if defined DEBUG_MAPPINGS_TO_STDERR || (defined DEBUG_TARGET && defined  DEBUG_PRINT_PBWT_MAPPINGS)
    altl = false; altr = false;
#endif
    }

    void set(const PBWTInterval &s) {
        i = s.i;
        j = s.j;
//        hfrom = s.hfrom;
//        hto = s.hto;
#if defined DEBUG_MAPPINGS_TO_STDERR || (defined DEBUG_TARGET && defined  DEBUG_PRINT_PBWT_MAPPINGS)
    altl = s.altl; altr = s.altr;
#endif
    }

    inline void setEmpty() {
        i = 0;
        j = -1;
#if defined DEBUG_MAPPINGS_TO_STDERR || (defined DEBUG_TARGET && defined  DEBUG_PRINT_PBWT_MAPPINGS)
    altl = false; altr = false;
#endif
    }

    int getLeft() const { return i; }
    int getRight() const { return j; }
//    Haplotype getHapFrom() const { return hfrom; }
//    Haplotype getHapTo() const { return hto; }

    bool operator==(const PBWTInterval &other) const {
        return i == other.i && j == other.j;
    }

    bool operator>(const PBWTInterval &other) const {
        // quite lazy:
        // if one interval is empty, the range is at least negative and the other positive -> correct;
        // if both are empty, both ranges are negative and which one is larger does not matter anyway
        return j - i > other.j - other.i;
    }

private:
    // borders
    int i,j;
    // haplotype where the corresponding path was originally mapped from (path begin)
//    Haplotype hfrom;
    // haplotype where the interval was mapped to (path end)
//    Haplotype hto;

#if defined DEBUG_MAPPINGS_TO_STDERR || (defined DEBUG_TARGET && defined  DEBUG_PRINT_PBWT_MAPPINGS)
public:
    bool altl = false, altr = false;
#endif
};

#endif /* PBWTINTERVAL_H_ */
