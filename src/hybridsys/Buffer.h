/*
 *    Copyright (C) 2018-2021 by Lars Wienbrandt and Jan Christian Kässens,
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

#ifndef BUFFER_H
#define BUFFER_H

#include <vector>
#include <list>
#include <map>
#include <mutex>
#include <memory>

#ifdef USE_AD_FPGA
#include <admxrc3.h>
#include <adb3.h>
#endif

#include "BufferAllocator.h"

namespace hybridsys {

template<typename T>
class BufferFactory;

template<typename T, template <typename> class Alloc>
class Buffer
{
    friend BufferFactory<Buffer<T,Alloc>>;

public:
    using value_type = T;
    using allocator = Alloc<T>;

#ifdef USE_AD_FPGA
    using BufferHandle = ADMXRC3_BUFFER_HANDLE;
    static constexpr BufferHandle BufferHandleInvalidValue = (unsigned)ADB3_HANDLE_INVALID_VALUE;
#else
    using BufferHandle = void*;
    static constexpr BufferHandle BufferHandleInvalidValue = nullptr;
#endif

    explicit Buffer(size_t buffer_size_, BufferHandle handle_ = BufferHandleInvalidValue)
        : buffer_size(buffer_size_), handle(handle_), content_length(0)
    {
        buffer.resize(buffer_size_);
    }

    ~Buffer() {}

    T *getData() {
        return buffer.data();
    }

    const T *getData() const {
        return buffer.data();
    }

    void   setContentLength(size_t length) { content_length = length; }
    size_t getContentLength() const { return content_length; }

    size_t getSize() const {
        return buffer_size;
    }

    BufferHandle getHandle() { return handle; }

private:

    const size_t buffer_size;
    BufferHandle handle;
    std::vector<T, Alloc<T>> buffer;
    size_t content_length;

};

using FPGABuffer = Buffer<char, PageAlignedAllocator>;
using CUDABuffer = Buffer<char, CUDAAllocator>;

}


#endif // BUFFER_H
