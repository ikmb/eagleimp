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

#ifdef USE_AD_FPGA
#include <admxrc3.h>
#endif

#ifdef USE_CUDA_GPU
#include <cuda_runtime.h>
#endif

#include "DeviceCategory.h"

namespace hybridsys {


const char *FPGACategory::name() const noexcept {
    return "FPGA";
}

const char *GPUCategory::name() const noexcept {
    return "GPU";
}

std::string FPGACategory::message(int code) const {
#ifdef USE_AD_FPGA
    return std::string(ADMXRC3_GetStatusString(static_cast<ADMXRC3_STATUS>(code), false));
#else
    (void) code;
    return std::string("Undefined");
#endif
}
std::string GPUCategory::message(int code) const {
#ifdef USE_CUDA_GPU
    return std::string(cudaGetErrorString(static_cast<cudaError_t>(code)));
#else
    (void) code;
    return std::string("Undefined");
#endif
}

}
