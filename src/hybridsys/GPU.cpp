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

#include <iomanip>
#include <sstream>

#ifdef USE_CUDA_GPU
#include <cuda_runtime.h>
#endif

#include "GPU.h"

namespace hybridsys {

GPU::GPU(int bus_, int slot_)
    : bus(bus_), slot(slot_), serial(""), cudaIndex(-1)
{
#ifdef USE_CUDA_GPU
    // construct an NVML bus id, assuming function and domain to be zero
    std::stringstream bus_id;
    bus_id << "0000:" << std::hex << std::setfill('0') << std::setw(2) << bus_ << ":" << slot_ << ".0";

    nvmlReturn_t ret = nvmlDeviceGetHandleByPciBusId(bus_id.str().c_str(), &nvmlHandle);
    if(ret != NVML_SUCCESS) {
        throw std::runtime_error(nvmlErrorString(ret));
    }

    // query for the serial number
    std::array<char, NVML_DEVICE_SERIAL_BUFFER_SIZE> serial_s;
    nvmlDeviceGetSerial(nvmlHandle, serial_s.data(), serial_s.size());
    this->serial = std::string(serial_s.data());

    // map a CUDA (non-NVML) device index to a physical ID
    int numGPUs = 0;
    cudaGetDeviceCount(&numGPUs);
    for(int i=0; i < numGPUs; i++) {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);
        if(props.pciBusID == bus_ && props.pciDeviceID == slot_)
            cudaIndex = i;
    }
    /* It is expected that NVML already bails if there is no CUDA device
       with that bus/device specification, so no additional checks here. */
#else
    throw std::runtime_error("CUDA libraries not available!");
#endif
}

int GPU::getBus() { return bus; }

int GPU::getSlot() { return slot; }

const std::string &GPU::getSerialNumber() { return serial; }

int GPU::getIndex() const { return cudaIndex; }

}
