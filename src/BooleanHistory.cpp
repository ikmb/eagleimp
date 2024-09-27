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

#include <vector>
#include <utility>
#include <iostream>

#include <stdint.h>

#include "StatusFile.h"

#include "BooleanHistory.h"

BooleanHistory::BooleanHistory(size_t historyLength) : histlast(historyLength-1) {
    if (historyLength == 0 || historyLength > DATA_VECTOR_SIZE * DATA_TYPE_BITS) {
        StatusFile::addError("BooleanHistory(len) called with an invalid history length.");
        exit(EXIT_FAILURE);
    }
    size_t last_word_idx = (historyLength-1)/DATA_TYPE_BITS;
    size_t mask_last_word;
    if (historyLength % DATA_TYPE_BITS == 0)
        mask_last_word = 0ull;
    else
        mask_last_word = ~((1ull << (historyLength%DATA_TYPE_BITS))-1ull);
    for (size_t i = 0; i < DATA_VECTOR_SIZE; i++) {
        if (i < last_word_idx)
            fillmask[i] = 0ull;
        else if(i == last_word_idx)
            fillmask[i] = mask_last_word;
        else
            fillmask[i] = ~0ull;
    }

    clear();
}

void BooleanHistory::push_new(bool b) { // like "push_front()", because raw comparisons rely on the most recent entries being stored at data word 0 from bit 0
    uint64_t bit = b ? 1ull : 0ull;

    // need to left shift the whole vector
    uint64_t rotbit = bit; // data will be inserted at lowest position
    for (uint64_t &d : data) { // bet this is fastest since compiler can optimize here
        uint64_t rotbitnew = d >> (DATA_TYPE_BITS-1); // the two most significant bits
        d = (d << 1) | rotbit;
        rotbit = rotbitnew;
    }
    if (last != histlast)
        last++;

    // ensure unused bits to be in a defined state
    for (size_t i = 0; i < DATA_VECTOR_SIZE; i++)
        data[i] |= fillmask[i];
}


bool BooleanHistory::at(size_t pos) const {
    size_t index = pos / DATA_TYPE_BITS;
    size_t bitindex = pos % DATA_TYPE_BITS;
    uint64_t mask = 1ull << bitindex;
    bool b = false;
    if (data[index] & mask)
        b = true;
    return b;
}

void BooleanHistory::clear() {
    last = -1;
    data.fill(0);

    // ensure unused bits to be in a defined state
    for (size_t i = 0; i < DATA_VECTOR_SIZE; i++)
        data[i] |= fillmask[i];
}

bool BooleanHistory::isEqual(const BooleanHistory &other) const {
    // due to a defined fill mask, comparison of all data components is sufficient
    bool ret = true;
    for (size_t i = 0; i <= DATA_VECTOR_SIZE; i++) // I believe this is fastest since the compiler can optimize here
        ret = ret && data[i] == other.data[i];
    return ret;
}

// returns 2*DATA_TYPE_BITS if the histories are completely equal, even if the history size is smaller!
size_t BooleanHistory::numEqual(const BooleanHistory &other) const {
    uint64_t t = data[0] ^ other.data[0];
    if (t != 0ull)
        return __builtin_ctzll(t);
//    return DATA_TYPE_BITS;
    else { // only if more than DATA_TYPE_BITS sites are equal
        t = data[1] ^ other.data[1];
        if (t != 0ull)
            return DATA_TYPE_BITS + __builtin_ctzll(t);
        else
            return 2*DATA_TYPE_BITS;
    }
}

//// DEBUG
//void BooleanHistory::printPath(ostream &out, int highlight1, int highlight2) const {
//    // implemented very slow, but used for debug purposes only anyway
//    for (int i = last; i >= 0; i--) {
//        if (i == highlight1 || i == highlight2)
//            out << " >";
//        out << at(i);
//        if (i == highlight1 || i == highlight2)
//            out << "< ";
//    }
//}
//// __DEBUG
//
//// DEBUG
//int64_t BooleanHistory::dbgCompare(const BooleanHistory &other, int delta) const {
//    uint64_t bitmask = (1ull << delta) - 1ull; // delta '1's, delta mustn't be greater than 63!
//    uint64_t bits = data[0] & bitmask;
//    uint64_t otherbits = other.data[0] & bitmask;
//    return (int64_t)dbgBitReverse(bits, delta) - (int64_t)dbgBitReverse(otherbits, delta);
//}
//// __DEBUG
//
//// DEBUG
//uint64_t dbgBitReverse(uint64_t data, size_t delta) {
//    uint64_t d = data;
//    uint64_t ret = 0ull;
//    for (size_t i = 0; i < delta; i++) {
//        ret <<= 1;
//        ret |= d & 1ull;
//        d >>= 1;
//    }
//    //cout << hex << data << "->" << ret << dec << endl;
//    return ret;
//}
//// __DEBUG


