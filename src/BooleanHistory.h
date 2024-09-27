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

#ifndef BOOLEANHISTORY_H_
#define BOOLEANHISTORY_H_

#include <iostream>
#include <utility>
#include <vector>
#include <array>
#include <cassert>

using namespace std;

//// DEBUG
//// slow bit-wise reverse
//uint64_t dbgBitReverse(uint64_t data, size_t delta);
//// __DEBUG

class BooleanHistory {
public:
    using data_type = uint64_t;
    static const size_t DATA_TYPE_BITS = sizeof(data_type) * 8;
    static const size_t DATA_VECTOR_SIZE = 2;

    BooleanHistory(size_t historyLength);
    BooleanHistory(const BooleanHistory&) = default;
    BooleanHistory(BooleanHistory&&) = default;
    BooleanHistory &operator=(const BooleanHistory&) = default;
    BooleanHistory &operator=(BooleanHistory&&) = default;

    void push_new(bool b); // pushes a new boolean
    bool at(size_t pos) const; /* no length check for speed reasons!! */
    bool operator[](size_t pos) const { return at(pos); }
    bool newest() const { return at(0); }
    bool oldest() const { return at(last); }
    size_t size() const { return static_cast<size_t>(last+1); }
    void clear();

    // returns true if both histories are completely equal
    bool isEqual(const BooleanHistory &other) const;

    // compares up to the last 128 sites of the path and returns min(128, the number of equal recent data slots (which can be more than the history size if all sites are equal!)).
    size_t numEqual(const BooleanHistory &other) const;

//    // DEBUG
//    void printPath(ostream &out, int highlight1, int highlight2 = -1) const;
//    // __DEBUG

//    // DEBUG
//    // for order check of the most recent delta sites
//    // returns:
//    // negative if this history < the other history
//    // 0 if equal
//    // positive if this history > the other history
//    int64_t dbgCompare(const BooleanHistory &other, int delta) const;
//    // __DEBUG

private:

    inline size_t numEqualInt(uint64_t a, uint64_t b, uint64_t mask) const;

    array<data_type, DATA_VECTOR_SIZE> fillmask; // is used to fill the unused bits in the data vector, for faster comparison
    array<data_type, DATA_VECTOR_SIZE> data;

    int last;
    int histlast;
};

#endif /* BOOLEANHISTORY_H_ */
