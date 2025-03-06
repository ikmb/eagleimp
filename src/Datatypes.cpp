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

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "Datatypes.h"

void GenotypeVector::push_back(const Genotype &g) {
    if (len % (GVUNITWORDS*DATATYPEBITS) == 0) {
        data_type is0 = (g == Genotype::HomRef || g == Genotype::Miss) ? 1ull : 0ull;
        data_type is2 = (g == Genotype::HomAlt || g == Genotype::Miss) ? 1ull : 0ull;
        data.push_back(is0);
        data.push_back(is2);

        // extend with a complete unit
        for (size_t i = 1; i < GVUNITWORDS; i++) {
            data.push_back(0ull); // is0
            data.push_back(0ull); // is2
        }

        nextpushword_is0 = data.data() + 2*(len/DATATYPEBITS);
        nextpushword_is2 = data.data() + 2*(len/DATATYPEBITS) + 1;
        nextpushbit = 2ull; // had to be 1ull before

    } else {
        data_type bit0 = (g == Genotype::HomRef || g == Genotype::Miss) ? nextpushbit : 0ull;
        data_type bit2 = (g == Genotype::HomAlt || g == Genotype::Miss) ? nextpushbit : 0ull;

        *nextpushword_is0 |= bit0;
        *nextpushword_is2 |= bit2;
        if(nextpushbit == HIGHESTPUSHBIT) {
            nextpushbit = 1ull;
            nextpushword_is0 += 2;
            nextpushword_is2 += 2;
        } else
            nextpushbit <<= 1;
    }
    len++;
}

void GenotypeVector::set(size_t pos, const Genotype &g) {
    size_t index = pos / DATATYPEBITS;
    size_t bitindex = pos % DATATYPEBITS;
    data_type bit = 1ull << bitindex;
    switch (g) {
    case Genotype::HomRef:
        data[2*index]   |= bit;  // set "is 0"
        data[2*index+1] &= ~bit; // clear "is 2"
        break;
    case Genotype::Het:
        data[2*index]   &= ~bit;
        data[2*index+1] &= ~bit;
        break;
    case Genotype::HomAlt:
        data[2*index]   &= ~bit;
        data[2*index+1] |= bit;
        break;
    case Genotype::Miss:
        data[2*index]   |= bit;
        data[2*index+1] |= bit;
        break;
    }
}

Genotype GenotypeVector::at(size_t pos) const {
    size_t index = pos / DATATYPEBITS;
    size_t bitindex = pos % DATATYPEBITS;
    data_type mask = 1ull << bitindex;
    if ((mask & data[2*index]) & data[2*index+1])
        return Genotype::Miss;
    if (mask & data[2*index])
        return Genotype::HomRef;
    if (mask & data[2*index+1])
        return Genotype::HomAlt;
    return Genotype::Het;
}

void GenotypeVector::reserve(size_t size) {
    size_t datasize = GVUNITWORDS * (size / (GVUNITWORDS*DATATYPEBITS) + (size % (GVUNITWORDS*DATATYPEBITS) ? 1 : 0));
    data.reserve(2*datasize); // 2 words: is0 and is2
}

void GenotypeVector::clear() {
    data.clear();
    len = 0;
}

size_t GenotypeVector::compare_withPreInit(const BooleanVector &hapvec) const {
    size_t nerr = 0;
    const auto &hdata = hapvec.getData();
    size_t limit = std::min(data.size()/2, hapvec.getDataSize()/sizeof(BooleanVector::data_type));
    for (size_t i = 0; i < limit; i++) {
        data_type tmp = (data[2*i] & hdata[i]) | (data[2*i+1] & ~hdata[i]);
        nerr += __builtin_popcountll(tmp);
    }
    return nerr;
}

void GenotypeVector::getInconsistencies(BooleanVector &inconret, const BooleanVector &hapvec) const {
    const data_type *hdata = hapvec.getData();
    data_type* retdata = inconret.getData();

    size_t limit = std::min(data.size()/2, hapvec.getDataSize()/sizeof(BooleanVector::data_type));
    for (size_t i = 0; i < limit; i++) {
        data_type tmp = (data[2*i] & hdata[i]) | (data[2*i+1] & ~hdata[i]);
        retdata[i] = tmp;
    }
}

// resizes the vector to keep only the last n elements
// also takes care of the memory to be erased properly (i.e. now unused areas are set to zero again)
void GenotypeVector::keepLast(size_t n) {
    // DEBUG
    if (n > len)
        cerr << "WARNING: GenotypeVector: tried to keep more elements (" << n << ") than the current size (" << len << ")." << endl;

    // keeping nothing is the same as clearing the vector
    if (n == 0) {
        clear();
        return;
    }

    // only need to proceed if n is lower than the number of stored elements
    if (n < len) {
        size_t lastidx = (len-1) / DATATYPEBITS;
        size_t lastbitidx = (len-1) % DATATYPEBITS;
        size_t new_lastidx = (n-1) / DATATYPEBITS;
        size_t new_lastbitidx = (n-1) % DATATYPEBITS;
        int64_t worddiff = lastidx - new_lastidx; // will be positive or zero

        // shift data words to align kept data correctly
        if (new_lastbitidx < lastbitidx) { // need a right shift
            size_t shift = lastbitidx - new_lastbitidx;
            data_type mask = (1ull << shift)-1; // to mask the carry-over bits
            data_type carry_is0 = 0; // for the bits to shift in
            data_type carry_is2 = 0; // for the bits to shift in
            for (int64_t i = lastidx; i >= worddiff; i--) {
                // store the bits that will be shifted out
                data_type new_carry_is0 = data[2*i] & mask;
                data_type new_carry_is2 = data[2*i+1] & mask;
                data[2*i] = (data[2*i] >> shift) | (carry_is0 << (DATATYPEBITS-shift));
                data[2*i+1] = (data[2*i+1] >> shift) | (carry_is2 << (DATATYPEBITS-shift));
                carry_is0 = new_carry_is0;
                carry_is2 = new_carry_is2;
            }
        } else if (new_lastbitidx > lastbitidx) { // need a left shift
            size_t shift = new_lastbitidx - lastbitidx;
            data_type mask = ((1ull << shift)-1) << (DATATYPEBITS-shift); // to mask the carry-over bits
            data_type carry_is0 = data[2*(worddiff-1)] & mask; // directly take the carry bits for the first word (note that worddiff will be positive here)
            data_type carry_is2 = data[2*(worddiff-1)+1] & mask; // directly take the carry bits for the first word (note that worddiff will be positive here)
            for (size_t i = worddiff; i <= lastidx; i++) {
                // store the bits that will be shifted out
                data_type new_carry_is0 = data[2*i] & mask;
                data_type new_carry_is2 = data[2*i+1] & mask;
                data[2*i] = (data[2*i] << shift) | (carry_is0 >> (DATATYPEBITS-shift));
                data[2*i+1] = (data[2*i+1] << shift) | (carry_is2 >> (DATATYPEBITS-shift));
                carry_is0 = new_carry_is0;
                carry_is2 = new_carry_is2;
            }
        } // else: new_lastbitidx == lastbitidx: no shift required

        // copy kept data to the front and erase previously used area (only if necessary)
        if (worddiff) {
            // copy
            for (size_t d = 0, s = worddiff; d <= new_lastidx; d++,s++) {
                data[2*d] = data[2*s];
                data[2*d+1] = data[2*s+1];
            }
            // erase == resize vector + extend to a complete unit
            data.resize(2*(new_lastidx+1));
            for (size_t i = (new_lastidx%GVUNITWORDS)+1; i < GVUNITWORDS; i++) {
                data.push_back(0ull); // is0
                data.push_back(0ull); // is2
            }
        }

        // update
        len = n;
        nextpushbit = 1ull << (n % DATATYPEBITS);
        nextpushword_is0 = data.data() + 2*(len/DATATYPEBITS);
        nextpushword_is2 = data.data() + 2*(len/DATATYPEBITS) + 1;
    }
}


ostream &operator<<(ostream &o, const Haplotype &hap) {
    o << (hap == Haplotype::Ref ? 0 : 1);
    return o;
}

ostream &operator<<(ostream &o, const Genotype &gen) {
    switch (gen) {
    case Genotype::HomRef:
        o << "0";
        break;
    case Genotype::Het:
        o << "1";
        break;
    case Genotype::HomAlt:
        o << "2";
        break;
    default: // Miss
        o << "X";
    }
    return o;
}

ostream &operator<<(ostream &o, const Constraint &c) {
    switch (c.type) {
    case Constraint::Type::HomRef:
        o << "Hom0";
        break;
    case Constraint::Type::HomAlt:
        o << "Hom1";
        break;
    case Constraint::Type::HetKeep:
        o << "HetK";
        break;
    case Constraint::Type::HetSwitch:
        o << "HetS";
        break;
    default: // no constraint
        o << "NoC";
        break;
    }
    o << "(" << c.refidx << ")";
    return o;
}

