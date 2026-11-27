/*
 *    Copyright (C) 2018-2025 by Lars Wienbrandt,
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

#include "RingBuffer.h"

bool RingBufferBool::at(size_t pos) const {
    if (pos >= len)
        cerr << "WARNING: RingBuffer: Accessing element out of range! size: " << len << " idx: " << pos << endl;

    if (pos+start < capacity)
        return buffer[pos+start];
    else
        return buffer[pos+start-capacity];
}

RingBufferBool::iterator RingBufferBool::begin() { return RingBufferBoolIterator(this, 0); }
RingBufferBool::iterator RingBufferBool::end() { return RingBufferBoolIterator(this, len); }
