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

#ifndef GPUKERNELS_H_
#define GPUKERNELS_H_

#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>

#include "GPUEngine.h"

#define checkCUDAError(err) { checkCUDAError_unwrapped((err), __FILE__, __LINE__); }
static void checkCUDAError_unwrapped(cudaError_t code, const char *file, int line) {
    if(code != cudaSuccess) {
        std::stringstream ss;
        ss << file << "(" << line << ")";
        std::string file_and_line;
        ss >> file_and_line;
        throw thrust::system_error(code, thrust::cuda_category(), file_and_line);
    }
}

// host wrapper for copying the constant symbols
__host__ void copyConstantsToDevice(
        unsigned K,
        unsigned Kwords,
        unsigned M,
        unsigned NrefhapsAD,
        size_t NrefhapsADwords,
        size_t tgtGTSize,
        size_t tgtSplitSize,
        size_t tgtBHSize);

// add inconsistencies to destination vector dependent on current genotype, existing inconsistencied will not be removed, size is taken from devNrefhapsADwords
__device__ inline void mergeInconsistencies(data_type* dest, const data_type* src, bool gtis0, bool gtis2, unsigned words);

// as memcpy with devNrefhapsADwords, but determining count0 on-the-fly
__device__ inline void copyAll(data_type* dest, const data_type* src, int* count0);

// the destination vector will consist only of the data found in src at the positions where besthaps is flagged with 1, size is taken from devNrefhapsAD
// zero count will be determined on-the-fly and stored at *count0
__device__ inline void copyBestHaps(data_type* dest, const data_type* src, const data_type* besthaps, int* count0);

// kernel creating the condensed reference for each target
__global__ void crefTKernel(
        const data_type* devRefTData,
        const uint32_t* devCrefTOffsets,
        const data_type* devTgtData,
        const data_type* devBHData,
        data_type* devCrefTData,
        data_type* devIncFullData,
        int* devCount0Data,
        int ntargets);

// kernel creating the compressed PBWT for each target:
// if there are more than 1 threads in a block, the creation is parallelized over the sites
// using an overlap. Result may thus not reflect the correct PBWT but a good approximation.
__global__ void pbwtKernel(
        const uint32_t* devCrefTOffsets,
        const data_type* devCrefTData,
        const int* devCount0Data,
        uint32_t* devPBWTData,
        uint32_t* devPermData,
        int ntargets,
        bool reverse);

// kernel creating the compressed PBWT for each target with parallelization over K:
// Parts are divided over K, the number of parts reflects the number of threads in a block.
// The results are merged afterwards and the offset fields correctly filled.
// Use a power of 2 for the number of threads per block!
__global__ void pbwtPartedKernel(
        const uint32_t* devCrefTOffsets,
        const data_type* devCrefTData,
        uint32_t* devPBWTData,
        uint32_t* devPermData,
        int ntargets,
        bool reverse);

__device__ inline void buildComprPBWTPart(
        const data_type* crefTData,
        uint32_t* comprPBWT,
        uint32_t* permData,
        int Mpbwt,
        int kpart,
        bool reverse);

// ptr is incremented by 2 if mask gets an overflow
__device__ inline void rdincrement(uint32_t* &ptr, uint32_t &mask, uint32_t &ptrval);

// ptr0 is incremented by 2 and ptr1 by 1 if mask gets an overflow
__device__ inline void rdincrement(uint32_t* &ptr0, uint32_t* &ptr1, uint32_t &mask, uint32_t &ptrval0, uint32_t &ptrval1);

// ptr0 is incremented by 2 and ptr1 by 1 if mask gets an overflow
__device__ inline void wrincrement(uint32_t* &ptr0, uint32_t* &ptr1, uint32_t &mask, uint32_t &ptrval0, uint32_t &ptrval1);

// parts to merge have to be located next to each other
// koff: beginning of parts to be merged (multiple of 32!)
// k0: length (in bits) of part 0 (multiple of 32!)
// k1: length (in bits) of part 1 (multiple of 32!)
//     NOTE: for the last part koff + k0 + k1 may exceed K, so be sure to have included the correct padding!
__device__ inline void merge(uint32_t* comprPBWT, uint32_t *orgspace, int Mpbwt, int koff, int k0, int k1, int pipedepth, int totaldepth);

// generate the offset information in the merged cPBWT. if 'odd' is set, the reference information is stored in the offset field, so this function will copy
// it to the data field while calculating the offset
__device__ inline void fillOffsets(uint32_t* comprPBWT, int Mpbwt, int pipedepth);


#endif /* GPUKERNELS_H_ */
