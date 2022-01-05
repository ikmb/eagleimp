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

#ifndef HAPLOTYPEHISTORY_H_
#define HAPLOTYPEHISTORY_H_

#include <iostream>
#include <utility>

#include "Datatypes.h"
#include "BooleanHistory.h"

class HaplotypeHistory {
public:
    HaplotypeHistory(size_t historyLength) : data(historyLength) {}
    HaplotypeHistory(const HaplotypeHistory&) = default;
    HaplotypeHistory(HaplotypeHistory&&) = default;
    HaplotypeHistory &operator=(const HaplotypeHistory&) = default;
    HaplotypeHistory &operator=(HaplotypeHistory&&) = default;

    // pushes a new haplotype pair (mat,pat)
    void push_new(const Haplotype &hap) { data.push_new(hap == Haplotype::Alt); }
    /* no length check for speed reasons!! */
    Haplotype at(size_t pos) const { return convertBool2Hap(data.at(pos)); }
    Haplotype operator[](size_t pos) const { return at(pos); }
    Haplotype newest() const { return convertBool2Hap(data.newest()); }
    Haplotype oldest() const { return convertBool2Hap(data.oldest()); }
    size_t size() const { return data.size(); }
    void clear() { data.clear(); }

    // returns true if the complete history is equal
    bool isEqual(const HaplotypeHistory &other) const { return data.isEqual(other.data); }

    // compares up to the last 128 sites of the path and returns min(128, the number of equal recent data slots (which can be more than the history size if all sites are equal!)).
    size_t numEqual(const HaplotypeHistory &other) const { return data.numEqual(other.data); }

//    // DEBUG
//    void printPath(ostream &out, int highlight1, int highlight2 = -1) const { data.printPath(out, highlight1, highlight2); }
//    // __DEBUG
//
//    // DEBUG
//    // for order check of the most recent delta sites
//    // returns:
//    // negative if this history < the other history
//    // 0 if equal
//    // positive if this history > the other history
//    int64_t dbgCompare(const HaplotypeHistory &other, int delta) const { return data.dbgCompare(other.data, delta); }
//    // __DEBUG

private:

    inline Haplotype convertBool2Hap(bool b) const { return b ? Haplotype::Alt : Haplotype::Ref; }

    BooleanHistory data;
};

#endif /* HAPLOTYPEHISTORY_H_ */
