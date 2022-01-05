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

#ifndef GPU_H
#define GPU_H

#ifdef USE_CUDA_GPU
extern "C" {
#include <nvml.h>
}
#endif

#include "Device.h"

namespace hybridsys {

class GPU : Device
{
public:
    GPU(int bus, int slot);

    // Device interface
    int getBus() override;
    int getSlot() override;
    const std::string& getSerialNumber() override;

    // Convert to device indices to be used with CUDA functions
    int getIndex() const override;

private:
    int bus;
    int slot;
    std::string serial;

    int cudaIndex;
#ifdef USE_CUDA_GPU
    nvmlDevice_t nvmlHandle;
#endif
};

}

#endif // GPU_H
