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

#ifndef GPUENGINE_H_
#define GPUENGINE_H_

#ifdef USE_CUDA_GPU

#include <cstdint>
#include <vector>
#include <sstream>

#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>
//extern "C" {
#include <cuda_runtime.h>
//#include <cuda.h>
//}

using namespace std;

using data_type = uint64_t;

class GPUEngine
{
public:

    GPUEngine() = delete;
    GPUEngine(int index);
    ~GPUEngine(){};

    // initializes the GPU, i.e. uploading required constants and reference to device, reserve space for src data and result
    void initialize(size_t K, size_t Kwords, size_t M, size_t Nrefhapscorr, size_t NrefhapsAD, size_t refTCapacity, const data_type* refTData);

    // returns the number of targets per block
    size_t getBlockSize() { return blocksize; }

    void runCrefPBWTKernel(
            const vector<const data_type*> &tgtData,
            const vector<const data_type*> &splitdata,
            const vector<const data_type*> &besthapdata,
            const vector<uint32_t> &crefTOffsets,
            vector<uint32_t*> &cPBWT,
            vector<vector<int>*> &count0s);

    // free space on device
    void free();

    int getCUDAIndex() const { return index; }

private:

    // global GPU index
    int index;

    // GPU specs
    cudaDeviceProp gpuProps;

    // the following will be set when initializing the engine
    size_t K;
    size_t Kwords; // (K+padding)/sizeof(data_type)
    size_t M;
    size_t Nrefhapscorr; // the number of haps that will be analyzed by each target
    size_t NrefhapsAD; // all haps in reference data block (including diploid encoding of haploids, depending on iteration with or without ALL phased targets)
    size_t blocksize; // how many targets per block?
    size_t crefThreadsPerTgt;
    size_t pbwtThreadsPerTgt;

    size_t tgtGTSize = 0;
    size_t tgtSplitSize = 0;
    size_t tgtBestHapSize = 0;

    // device addresses, space will be reserved when initializing the engine
    data_type* devRefTData;
    uint32_t*  devCrefTOffsets;

    data_type* devTgtData;
    data_type* devBHData;
    data_type* devTgtTmpData;

};

#endif // USE_CUDA_GPU

#endif /* GPUENGINE_H_ */
