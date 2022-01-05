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

        nextpushword_is0 = &(data.data()[2*(len/DATATYPEBITS)]);
        nextpushword_is2 = &(data.data()[2*(len/DATATYPEBITS)+1]);
        nextpushbit = 2ull; // had to be 1ull before

    } else {
        data_type bit0 = (g == Genotype::HomRef || g == Genotype::Miss) ? nextpushbit : 0ull;
        data_type bit2 = (g == Genotype::HomAlt || g == Genotype::Miss) ? nextpushbit : 0ull;
//        // is0
//        data[2*(len/DATATYPEBITS)] |= bit0 << (len % DATATYPEBITS); // ok if bit is zero
//        // is2
//        data[2*(len/DATATYPEBITS)+1] |= bit2 << (len % DATATYPEBITS); // because all is pre-initialized with zeros
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

