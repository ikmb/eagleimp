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

#ifndef RINGBUFFER_H_
#define RINGBUFFER_H_

#include <cstdint>
#include <cstring>
#include <vector>
#include <iostream>

using namespace std;

template<typename T>
class RingBufferIterator;

template<typename T>
class RingBuffer {

public:
    RingBuffer(){
        reset(16); // some initial capacity
    }

    RingBuffer(size_t _capacity, size_t _size = 0) {
        reset(_capacity, _size);
    }

    ~RingBuffer() = default;

    // resets the buffer with given capacity and optionally sets the initial size
    void reset(size_t _capacity, size_t _size = 0) {
        if (_size > _capacity) {
            capacity = _size;
        } else
            capacity = _capacity;

        buffer.clear();
        buffer.resize(capacity);

        len = _size;
        start = 0;
        last = _capacity == _size ? 0 : _size;
    }

    // insert an element at the end, increase (double) capacity if required
    void push_back(const T& e) {
        // increase capacity if required
        if (len == capacity) {
            buffer.resize(2*capacity); // we double here to ensure there's enough new space in case we need to copy the whole old vector now
            // copy all elements from 0 to start to the new allocated area at the end of the old vector
            // NOTE: at this point is last == start, and the old vector was full.
            for (size_t i=0; i < start; i++) {
                buffer[capacity+i] = buffer[i];
            }
            last = capacity + start;
            capacity *= 2;
        }
        // insert element
        buffer[last] = e;
        last++;
        if (last == capacity)
            last = 0;
        len++;
    }

    // as alternative to lvalue back(): modify the last element (no check for empty container!)
    void set_back(const T& e) {
        if (last > 0)
            buffer[last-1] = e;
        else
            buffer[capacity-1] = e;
    }

    // get a copy of an element, no check if this was already inserted. (faster!)
    T& at(size_t pos) {
        if (pos >= len)
            cerr << "WARNING: RingBuffer: Accessing element out of range! size: " << len << " idx: " << pos << endl;

        if (pos+start < capacity)
            return buffer[pos+start];
        else
            return buffer[pos+start-capacity];
    }

    // get a const reference to an element, no check if this was already inserted. (faster!)
    const T& at(size_t pos) const {
        if (pos >= len)
            cerr << "WARNING: RingBuffer: Accessing element out of range! size: " << len << " idx: " << pos << endl;

        if (pos+start < capacity)
            return buffer[pos+start];
        else
            return buffer[pos+start-capacity];
    }

    T& back() {
        if (len == 0)
            cerr << "WARNING: RingBuffer: Accessing last element of empty buffer!" << endl;
        return at(len-1);
    }

    const T& back() const {
        return at(len-1);
    }

    T& operator[](size_t pos) {
        return at(pos);
    }

    const T& operator[](size_t pos) const {
        return at(pos);
    }

    size_t size() const {
        return len;
    }

    void clear() {
        start = 0;
        last = 0;
        len = 0;
    }

    // reduces the buffer size to keep only the last n elements
    void keepLast(size_t n) {
        if (n > len)
            cerr << "WARNING: RingBuffer: Attempting to keep more elements than available." << endl;

        if (n < len) {
            start = last >= n ? last-n : capacity+last-n;
            len = n;
        }
    }

    // sets the buffer size to the given value, but doesn't touch the underlying data!
    // if increased, the end pointer is shifted, data is not touched, so it may be undefined!
    // if decreased, data is removed (but not destroyed) from the back by moving the end pointer.
    // no check if data is increased over capacity!
    void setSize(size_t n) {
        last = (start + n) >= capacity ? (start + n - capacity) : (start + n);
        len = n;
    }

    // subtracts (operator-) n from every element
    void sub(T n) {
        if (last > start) {
            for (size_t p = start; p < last; p++)
                buffer[p] -= n;
        } else {
            for (size_t p = 0; p < last; p++)
                buffer[p] -= n;
            for (size_t p = start; p < capacity; p++)
                buffer[p] -= n;
        }
    }

    typedef RingBufferIterator<T> iterator;
    iterator begin() { return RingBufferIterator<T>(this, 0); }
    iterator end() { return RingBufferIterator<T>(this, len); }

protected:
    vector<T> buffer;
    size_t capacity = 0;
    size_t start = 0; // inclusive
    size_t last = 0;  // exclusive
    size_t len = 0;
};

template<typename T>
class RingBufferIterator {

public:
    RingBufferIterator() = delete;
    RingBufferIterator(RingBuffer<T>* rb_, size_t it_ = 0) : rb(rb_), it(it_) {}
    ~RingBufferIterator(){}
    RingBufferIterator(const RingBufferIterator<T>&) = default;

    // no boundary checks!
    RingBufferIterator<T>& operator++() { it++; return (*this); }
    RingBufferIterator<T>& operator--() { it--; return (*this); }
    RingBufferIterator<T>& operator+=(size_t i) { it+=i; return (*this); }
    RingBufferIterator<T>& operator-=(size_t i) { it-=i; return (*this); }
    bool operator==(const RingBufferIterator<T>& other) const { return other.rb == rb && other.it == it; };
    bool operator!=(const RingBufferIterator<T>& other) const { return other.rb != rb || other.it != it; };
    T& operator*() { return rb->at(it); }
    const T& operator*() const { return rb->at(it); }
    T* operator->() { return &(rb->at(it)); }

protected:
    RingBuffer<T>* rb;
    size_t it;
};




class RingBufferBoolIterator;

class RingBufferBool : RingBuffer<bool> {
public:

    // resets the buffer with given capacity and optionally sets the initial size
    void reset(size_t _capacity, size_t _size = 0) {
        RingBuffer::reset(_capacity, _size);
    }

    // insert an element at the end, increase (double) capacity if required
    void push_back(bool e) {
        RingBuffer::push_back(e);
    }

    // as alternative to lvalue back(): modify the last element (no check for empty container!)
    void set_back(bool e) {
        RingBuffer::set_back(e);
    }

    bool at(size_t pos) const;

    bool operator[](size_t pos) const {
        return at(pos);
    }

    bool back() const {
        return at(len-1);
    }

    size_t size() const {
        return RingBuffer::size();
    }

    void clear() {
        RingBuffer::clear();
    }

    // reduces the buffer size to keep only the last n elements
    void keepLast(size_t n) {
        RingBuffer::keepLast(n);
    }

    typedef RingBufferBoolIterator iterator;
    iterator begin();
    iterator end();
};

class RingBufferBoolIterator {

public:
    RingBufferBoolIterator() = delete;
    RingBufferBoolIterator(RingBufferBool* rb_, size_t it_ = 0) : rb(rb_), it(it_) {}
    ~RingBufferBoolIterator(){}
    RingBufferBoolIterator(const RingBufferBoolIterator&) = default;

    // no boundary check!
    RingBufferBoolIterator& operator++() { it++; return (*this); }
    RingBufferBoolIterator& operator--() { it--; return (*this); }
    RingBufferBoolIterator& operator+=(size_t i) { it+=i; return (*this); }
    RingBufferBoolIterator& operator-=(size_t i) { it-=i; return (*this); }
    bool operator==(const RingBufferBoolIterator& other) const { return other.rb == rb && other.it == it; };
    bool operator!=(const RingBufferBoolIterator& other) const { return other.rb != rb || other.it != it; };
    bool operator*() { return rb->at(it); }

protected:
    RingBufferBool* rb;
    size_t it;
};


#endif /* RINGBUFFER_H_ */
