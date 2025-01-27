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

#ifndef DATATYPES_H_
#define DATATYPES_H_

#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>
#include <array>
#include <iostream>

#include "StatusFile.h"
#include "utils.h"

// all PBWT mappings to stderr (independent of DEBUG_TARGET!)
//#define DEBUG_MAPPINGS_TO_STDERR
// create a small profile of PBWT mappings
//#define DEBUG_PROFILE_PBWT_MAPPINGS
// a lot of debug output for each step for given target range
//#define DEBUG_TARGET
// no additional debug output but phasing is restricted to below range
//#define DEBUG_TARGET_SILENT
// only few debug output per target in given range incl. per-site probabilities
//#define DEBUG_TARGET_LIGHT
// target range (all including)
#define DEBUG_TARGET_START 0
#define DEBUG_TARGET_STOP  63
//#define DEBUG_TARGET_START 112
//#define DEBUG_TARGET_STOP  211
// additional print-outs if DEBUG_TARGET is defined
//#define DEBUG_PRINT_PBWT_MAPPINGS
//#define DEBUG_PRINT_PATHS_AFTER_EXT
//#define DEBUG_PRINT_PATHS_AFTER_PRUNE
//#define DEBUG_PRINT_PATHS_AFTER_EXT_APP
//#define DEBUG_PRINT_PATHS_AFTER_MERGE
//#define DEBUG_PRINT_SORTED_BY_PROB
//#define DEBUG_PRINT_PREC
//#define DEBUG_PRINT_BESTHAPS
//#define DEBUG_PRINT_LOCATIONS
//#define DEBUG_PRINT_INCONS
//#define DEBUG_PRINT_IBD
//#define DEBUG_PRINT_PROBS_SUMMARY

//#define STOPWATCH


using namespace std;

using fp_type = float;

const size_t UNITWORDS = 8; // needs to be a multiple of 2 for GenotypeVector!!!

enum class Haplotype : bool {
    Ref = false,
    Alt = true
};

inline Haplotype operator!(Haplotype h) { return h == Haplotype::Ref ? Haplotype::Alt : Haplotype::Ref; }
inline bool toBool(Haplotype h) { return h == Haplotype::Ref ? false : true; }
inline Haplotype toHap(bool b) { return b ? Haplotype::Alt : Haplotype::Ref; }

enum class Genotype {
    HomRef = 0,
    Het    = 1,
    HomAlt = 2,
    Miss   = 3
};

class Constraint {
public:
    enum class Type {
        HomRef,
        HetKeep,
        HetSwitch,
        HomAlt,
        NoConstr
    } type;
    int refidx; // the reference index where to look at in the path history before(!) extension (i.e. 0 for the newest hap)
    // default: unconstrained
    Constraint() : type(Type::NoConstr), refidx(0) {}
    Constraint(Type type_, int refidx_) : type(type_), refidx(refidx_) {}
};

class PathProbability {
public:
    fp_type p;
    int scaleexp;
    PathProbability() : p(0.0), scaleexp(0) {}
    PathProbability(fp_type p_, int scaleexp_) : p(p_), scaleexp(scaleexp_) {}
};

class BooleanVector {
public:
    using data_type = std::uint64_t;
    static const size_t DATATYPEBITS = sizeof(data_type) * 8;

    BooleanVector() : len(0), capacity(0), data(NULL), nextpushbit(1ull), nextpushword(NULL) {};

    // initialize vector with no contents, data space not initialized!
    BooleanVector(data_type* data_, size_t capacity_) {
        init(data_, capacity_, 0ull);
    }

    // initialize vector to have given size, content not predefined
    BooleanVector(data_type* data_, size_t capacity_, size_t size) {
        init(data_, capacity_, size);
    }

    // initialize vector to have given size, filled with initValue
    BooleanVector(data_type* data_, size_t capacity_, size_t size, bool initValue) {
        init(data_, capacity_, size);
        setInitV(initValue);
    }

    inline void setData(data_type* data_, size_t capacity_, size_t size) {
        init(data_, capacity_, size);
    }

    inline void setDataAndInit(data_type* data_, size_t capacity_, size_t size, bool initValue) {
        init(data_, capacity_, size);
        setInitV(initValue);
    }

    // keeps current data space and capacity but sets current length to zero and clears the data space
    inline void clear() {
        init(data, capacity, 0);
        setInitV(false);
    }

    // keeps the underlying data as is and sets the internal size to the provided value
    inline void setSize(size_t size_) {
        init(data, capacity, size_);
    }

    // only works if capacity is set large enough and memory is preinitialized with zero! otherwise you risk a segfault or undefined behaviour, no length check!!!
    inline void push_back_withPreInit(bool b) {
        if (b)
            *nextpushword |= nextpushbit;

        if (nextpushbit == HIGHESTPUSHBIT) {
            nextpushbit = 1ull;
            nextpushword++;
        } else
            nextpushbit <<= 1;
        len++; // this is ok if bit is zero because all bits are supposed to be preinitialized with zeros

    }

    /* no length check for speed reasons!! */
    inline void set(size_t pos, bool b) {
        size_t index = pos / DATATYPEBITS;
        size_t bitindex = pos % DATATYPEBITS;
        data_type bit = 1ull << bitindex;
        if (b) // 1
            data[index] |= bit; // set
        else // 0
            data[index] &= ~bit; // clear
    }

    /* no length check for speed reasons!! */
    inline void setWithPreInit(size_t pos, bool b) {
        size_t index = pos / DATATYPEBITS;
        size_t bitindex = pos % DATATYPEBITS;
        if (b) // 1
            data[index] |= 1ull << bitindex; // set
    }

    // sets two positions at once, expects pre-initialization with zeros (false)
    /* no length check for speed reasons!! */
    /* only set at even positions! pos%2 must be zero! */
    inline void setPairWithPreInit(size_t pos, bool b1, bool b2) {
        size_t index = pos / DATATYPEBITS;
        size_t bitindex = pos % DATATYPEBITS;
        if (b1)
            data[index] |= 1ull << bitindex; // set
        if (b2)
            data[index] |= 2ull << bitindex; // set
    }

    // all bits AFTER this position are cleared, but the length of the vector is not changed
    inline void clearAfter(size_t pos) {
        // first position after provided
        size_t index = (pos+1) / DATATYPEBITS;
        size_t bitindex = (pos+1) % DATATYPEBITS;
        data_type mask = (1ull << bitindex) - 1ull;
        data[index] &= mask; // sets all bits after provided position to zero in this field
        // set all remaining fields to zero
        index++;
        for (; index < capacity / sizeof(data_type); index++)
            data[index] = 0ull;
    }

    /* no length check for speed reasons!! */
    inline bool at(size_t pos) const {
        size_t index = pos / DATATYPEBITS;
        size_t bitindex = pos % DATATYPEBITS;
        data_type mask = 1ull << bitindex;
        return data[index] & mask;
    }

    inline size_t size() const { return len; }
    // returns the current size of the data in the vector in bytes, rounded to a full unit
    inline size_t getDataSize() const { return roundToMultiple(len, UNITWORDS * sizeof(data_type) * 8) / 8; }
    inline size_t getCapacity() const { return capacity; }

    inline data_type* getData() { return data; }
    const inline data_type* getData() const { return data; }

    // comparison only up to common length, returns the number of differences
    inline size_t compare(const BooleanVector &hapvec) const {
        size_t nerr = 0;
        const data_type *hdata = hapvec.getData();
        size_t limit = std::min(len, hapvec.size());
        for (size_t i = 0; i < limit/sizeof(data_type); i++) {
            data_type tmp = data[i] ^ hdata[i]; // xor: all different positions are set
            nerr += __builtin_popcountll(tmp);
        }
        if (limit % sizeof(data_type)) {
            data_type mask = (1ull << (limit % sizeof(data_type))) - 1ull;
            data_type tmp = mask & (data[limit/sizeof(data_type)] ^ hdata[limit/sizeof(data_type)]); // xor: all different positions are set
            nerr += __builtin_popcountll(tmp);
        }
        return nerr;
    }

    // returns the number of false bits in this vector,
    // works only correctly if memory was pre-initialized with zeros!!!
    inline size_t getFalseCount_withPreInit() const {
        // simply count all set bits with popcount (if a bit is set, it must be included in the data, since unused bits are unset)
        int cnt1 = 0;
        for (size_t i = 0; i < divideRounded(len, 8*sizeof(data_type)); i++)
            cnt1 += __builtin_popcountll(data[i]);
        // return the length of the vector minus the count of all set bits
        return len - (size_t) cnt1;
    }

    /* no length check for speed reasons!! */
    inline bool operator[](size_t pos) const { return at(pos); }

    // resizes the vector to keep only the last n elements
    // also takes care of the memory to be erased properly (i.e. now unused areas are set to zero again)
    inline void keepLast(size_t n) {
        // DEBUG
        if (n > len)
            cerr << "WARNING: BooleanVector: tried to keep more elements (" << n << ") than the current size (" << len << ")." << endl;

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
            size_t worddiff = lastidx - new_lastidx; // will be positive or zero

            // shift data words to align kept data correctly
            if (new_lastbitidx < lastbitidx) { // need a right shift
                size_t shift = lastbitidx - new_lastbitidx;
                data_type mask = (1ull << shift)-1; // to mask the carry-over bits
                data_type carry = 0; // for the bits to shift in
                for (size_t i = lastidx; i >= worddiff; i--) {
                    // store the bits that will be shifted out
                    data_type new_carry = data[i] & mask;
                    data[i] = (data[i] >> shift) | (carry << (DATATYPEBITS-shift));
                    carry = new_carry;
                }
            } else if (new_lastbitidx > lastbitidx) { // need a left shift
                size_t shift = new_lastbitidx - lastbitidx;
                data_type mask = ((1ull << shift)-1) << (DATATYPEBITS-shift); // to mask the carry-over bits
                data_type carry = data[worddiff-1] & mask; // directly take the carry bits for the first word (note that worddiff will be positive here)
                for (size_t i = worddiff; i <= lastidx; i++) {
                    // store the bits that will be shifted out
                    data_type new_carry = data[i] & mask;
                    data[i] = (data[i] << shift) | (carry >> (DATATYPEBITS-shift));
                    carry = new_carry;
                }
            } // else: new_lastbitidx == lastbitidx: no shift required

            // copy kept data to the front and erase previously used area (only if necessary)
            if (worddiff) {
                // copy
                for (size_t d = 0, s = worddiff; d <= new_lastidx; d++,s++)
                    data[d] = data[s];
                // erase
                for (size_t i = new_lastidx+1; i <= lastidx; i++)
                    data[i] = 0;
            }

            // update
            len = n;
            nextpushbit = 1ull << (n % DATATYPEBITS);
            nextpushword = data + (n / DATATYPEBITS);
        }
    }

    inline void dump(ostream &out) const {
        size_t ac = 0;
        for (size_t i = 0; i < len; i++) {
            out << ' ';
            if (at(i)) {
                out << '1';
                ac++;
            } else
                out << '0';
        }
        out << " ac: " << ac << endl;
    }


private:
    inline void init(data_type* data_, size_t capacity_, size_t size) {
        if (capacity_ % (UNITWORDS * sizeof(data_type))) {
            StatusFile::addError("BooleanVector: capacity must be a multiple of UNITWORDS*sizeof(data_type)!");
            exit(EXIT_FAILURE);
        }
        capacity = capacity_;
        data = data_;
        len = size;
        nextpushbit = 1ull << (size % (8*sizeof(data_type)));
        nextpushword = data + (size / (8*sizeof(data_type)));
    }

    inline void setInitV(bool initValue) {
        int initv = initValue ? -1 : 0;
        memset(data, initv, capacity);
        // we have to clear the leading unused bits of the newest UNITWORDS data words if we initialized with ones
        if (initValue && (len % (UNITWORDS*DATATYPEBITS))) {
            data_type mask = (1ull << (len % DATATYPEBITS)) - 1ull; // all used bits are one, the others are zero
            size_t mask_word_idx = (len % (UNITWORDS*DATATYPEBITS)) / DATATYPEBITS; // relative word in the unit where the mask has to be applied
            size_t datasize = capacity / sizeof(data_type);
            mask_word_idx = datasize - UNITWORDS + mask_word_idx; // absolute data index for mask application
            data[mask_word_idx] = mask; // replace the last word with the mask (works, since the current word is 11...11)
            for (size_t wd = mask_word_idx+1; wd < datasize; wd++)
                data[wd] = 0ull; // clear the remaining leading words
        }
    }

    size_t len;
    size_t capacity; // capacity in number of bytes
    data_type *data;
    data_type nextpushbit;
    data_type *nextpushword;
    static const data_type HIGHESTPUSHBIT = 1ull << (8*sizeof(data_type)-1);
};


class GenotypeVector {
public:
    using data_type = std::uint64_t;
    static const size_t DATATYPEBITS = sizeof(data_type) * 8;
    static const size_t GVUNITWORDS = UNITWORDS; // for genotype vectors UNITWORDS needs to be a multiple of 2, since it internally stores gt data in two words: is0 and is2

    GenotypeVector() : len(0), nextpushbit(1ull), nextpushword_is0(NULL), nextpushword_is2(NULL) {}
    void push_back(const Genotype &g);
    void set(size_t pos, const Genotype &g); /* no length check for speed reasons!! */
    Genotype at(size_t pos) const; /* no length check for speed reasons!! */
    Genotype operator[](size_t pos) const { return at(pos); }
    void reserve(size_t size);
    size_t size() const { return len; }
    void clear();

    // comparison only up to common length
    // returns the number of inconsistencies (missing positions are always registered as an inconsistency!)
    // expects the BooleanVector to be pre-initialized with zeros!
    size_t compare_withPreInit(const BooleanVector &hapvec) const;

    // returns a BooleanVector indicating all inconsistent sites in inconret. The vector must at least have the capacity of hapvec size.
    void getInconsistencies(BooleanVector &inconret, const BooleanVector &hapvec) const;

    // resizes the vector to keep only the last n elements
    // also takes care of the memory to be erased properly (i.e. now unused areas are set to zero again)
    void keepLast(size_t n);

    // returns the data in the following encoding: even words are "is0", odd words are "is2".
    const vector<data_type> &getData() const { return data; }

private:
    size_t len;
    // encoding: even words are "is0", odd words are "is2"
    vector<data_type> data;

    data_type nextpushbit;
    data_type *nextpushword_is0;
    data_type *nextpushword_is2;
    static const data_type HIGHESTPUSHBIT = 1ull << (8*sizeof(data_type)-1);
};

ostream &operator<<(ostream&, const Haplotype&);
ostream &operator<<(ostream&, const Genotype&);
ostream &operator<<(ostream&, const Constraint&);

#endif /* DATATYPES_H_ */
