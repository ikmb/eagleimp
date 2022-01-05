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

#ifndef HISTORY_H_
#define HISTORY_H_

#include <stdexcept>
#include <vector>

#include "StatusFile.h"

template<class T>
class History {

public:
    explicit History(size_t historyLength) : histlen(historyLength), start(0)
    {
        if (historyLength == 0) {
            StatusFile::addError("History(len) called with an invalid history length (0).");
            exit(EXIT_FAILURE);
        }
        log2datasize = 63 - __builtin_clzll(historyLength);
        if ((1ull << log2datasize) != historyLength)
            log2datasize++; // round up
        data.clear();
        data.resize(1ull << log2datasize); // always space with the size a power of 2
    }
    History(const History&) = default;
    History(History&&) = default;
    History &operator=(const History&) = default;
    History &operator=(History&&) = default;

    void push_new(const T &t) { // like push_front()
        start = modulo(start - 1);
        data[start] = t;
    }

    inline T &at(size_t pos) { /* no length check for speed reasons!! */
        return data[modulo(pos+start)];
    }

    const inline T &at(size_t pos) const { /* no length check for speed reasons!! */
        return data[modulo(pos+start)];
    }

    T &operator[](size_t pos) { return at(pos); }
    const T &operator[](size_t pos) const { return at(pos); }
    const T &newest() const { return data[start]; }
    void clear() {
        start = 0;
    }

private:
    // checks bounds of an unchecked array position and returns the correct array index (basically a modulo operation)
    inline size_t modulo(size_t pos) const {
        return pos % (1ull << log2datasize); // cool, compiler recognizes the mask here!
    }

    size_t histlen;
    size_t start; // including
    std::vector<T> data;
    size_t log2datasize;

};

#endif /* HISTORY_H_ */
