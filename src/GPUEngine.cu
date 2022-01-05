#ifdef USE_CUDA_GPU

#include <sstream>
#include <iostream>
#include <algorithm>

#include "Args.h"

#include "GPUEngine.h"
#include "GPUHandler.h"
#include "GPUKernels.h"

GPUEngine::GPUEngine(int index_)
    : index(index_)
{
    // the following will be set when initializing the engine
    K = 0;
    Kwords = 0;
    M = 0;
    Nrefhapscorr = 0;
    NrefhapsAD = 0;
    blocksize = 0;
    crefThreadsPerTgt = 0;
    pbwtThreadsPerTgt = 0;

    devRefTData     = NULL;
    devCrefTOffsets = NULL;

    devTgtData      = NULL;
    devBHData       = NULL;
    devTgtTmpData  = NULL;

//    devSrcBuffer = NULL;
//    devDstBuffer = NULL;
////    devRevBuffer = NULL;
////
////    devNumConsBuf = NULL;
////    devCnt0Buf = NULL;
////
////    devAccMvec = NULL;
////    devTargetOffsets = NULL;
}

void GPUEngine::initialize(size_t K_, size_t Kwords_, size_t M_, size_t Nrefhapscorr_, size_t NrefhapsAD_, size_t refTCapacity, const data_type* refTData) {
    cudaSetDevice(index);

    K = K_;
    Kwords = Kwords_;
    M = M_;
    Nrefhapscorr = Nrefhapscorr_;
    NrefhapsAD = NrefhapsAD_;
    size_t NrefhapsADwords = refTCapacity / sizeof(data_type);

    checkCUDAError(cudaGetDeviceProperties(&gpuProps, index))

    cout << "GPU Probs: " << index << endl;
    cout << "MPU count: " << gpuProps.multiProcessorCount << endl;
    cout << "Max threads per block: " << gpuProps.maxThreadsPerBlock << endl;
    cout << "Max threads per MPU: " << gpuProps.maxThreadsPerMultiProcessor << endl;
    cout << "Shared memory per block: " << gpuProps.sharedMemPerBlock << endl;
    cout << "Shared memory per MPU: " << gpuProps.sharedMemPerMultiprocessor << endl;
    cout << "Total global memory: " << gpuProps.totalGlobalMem << endl;
    cout << "Mem bus width: " << gpuProps.memoryBusWidth << endl;

    // TODO according to the capabilities: determine number of targets per block and block dimension
    blocksize = 256;
    crefThreadsPerTgt = 512;
    pbwtThreadsPerTgt = 16; // must be a power of 2!!!

    // reserve space on device for reference
    checkCUDAError(cudaMalloc(&devRefTData, M * refTCapacity))

    checkCUDAError(cudaMalloc(&devCrefTOffsets, (blocksize+1)*sizeof(uint32_t))) // offsets for the destination data

    // reserve space for target input
    tgtGTSize = roundToMultiple(M*2, 8*sizeof(data_type)*UNITWORDS) / 8; // tgt genotype data
    tgtSplitSize = roundToMultiple(M, 8*sizeof(data_type)*UNITWORDS) / 8;   // tgt split site information
    tgtBestHapSize = refTCapacity; // tgt best haps information
    // TODO Will it be much faster, if target data per block data is available in transposed form?
    size_t tgtsrcsize = tgtGTSize + tgtSplitSize;
    checkCUDAError(cudaMalloc(&devTgtData, blocksize * tgtsrcsize)) // single malloc for all GTs and splits
    checkCUDAError(cudaMalloc(&devBHData, blocksize * tgtBestHapSize)) // best haps data
    // temporary space for collecting incon information between split sites (for cref creation)
    // and global permutation arrays (for PBWT creation) - requires 2*K ints
    // and origins information (for PBWT merging) - requires 2*log_2(pbwtthreads)*K bits, which is always less than the size for the global perm arrays
    size_t tmpsize = max(crefThreadsPerTgt * refTCapacity, 2 * K * sizeof(int));
    checkCUDAError(cudaMalloc(&devTgtTmpData, blocksize * tmpsize))

    // space for target output cannot be reserved yet since we don't know the number of split sites per target yet

    // copy constant symbols to device
    copyConstantsToDevice(K, Kwords, M, NrefhapsAD, NrefhapsADwords, tgtGTSize, tgtSplitSize, tgtBestHapSize);

    // copy reference to device
    checkCUDAError(cudaMemcpy(devRefTData, refTData, M*refTCapacity, cudaMemcpyHostToDevice))
}

void GPUEngine::runCrefPBWTKernel(
        const vector<const data_type*> &tgtData,
        const vector<const data_type*> &splitdata,
        const vector<const data_type*> &besthapdata,
        const vector<uint32_t> &crefTOffsets,
        vector<uint32_t*> &cPBWT,
        vector<vector<int>*> &count0s) {
    cudaSetDevice(index);

    int threadsPerBlock = crefThreadsPerTgt;
    int blocksPerGrid = blocksize;

//    // DEBUG
//    cout << "GPU crefT kernel: " << endl;
//    cout << " Block size:        " << blocksize << endl;
//    cout << " Threads per block: " << threadsPerBlock << endl;
//    cout << " Blocks per grid:   " << blocksPerGrid << endl;

    // allocate memory for return data
    data_type* devCrefTData;
    int* devCount0Data;
    checkCUDAError(cudaMalloc(&devCrefTData, crefTOffsets.back()*Kwords*sizeof(data_type)))
    checkCUDAError(cudaMalloc(&devCount0Data, crefTOffsets.back()*sizeof(int)))

    // copy destination offsets to device
    checkCUDAError(cudaMemcpy(devCrefTOffsets, crefTOffsets.data(), crefTOffsets.size()*sizeof(uint32_t), cudaMemcpyHostToDevice))

    // copy target data to device
    size_t currtgtoff = 0;
    size_t currbhoff = 0;
    for (size_t t = 0; t < tgtData.size(); t++) {
//        cout << t << "/" << tgtData.size() << ":" << endl;
//        cout << "tgtGT:      @" << hex << (unsigned long)(devTgtData+curroff) << ": " << tgtData[t] << " size: " << dec << tgtGTSize << endl;
        checkCUDAError(cudaMemcpy(devTgtData+currtgtoff, tgtData[t], tgtGTSize, cudaMemcpyHostToDevice))
        currtgtoff += tgtGTSize/sizeof(data_type);
//        cout << "tgtSplits:  @" << hex << (unsigned long)(devTgtData+curroff) << ": " << splitdata[t] << " size: " << dec << tgtSplitSize << endl;
        checkCUDAError(cudaMemcpy(devTgtData+currtgtoff, splitdata[t], tgtSplitSize, cudaMemcpyHostToDevice))
        currtgtoff += tgtSplitSize/sizeof(data_type);
//        cout << "tgtBesthap: @" << hex << (unsigned long)(devTgtData+curroff) << ": " << besthapdata[t] << " size: " << dec << tgtBestHapSize << endl;
        checkCUDAError(cudaMemcpy(devBHData+currbhoff, besthapdata[t], tgtBestHapSize, cudaMemcpyHostToDevice))
        currbhoff += tgtBestHapSize/sizeof(data_type);
    }

    // run kernel
    cudaThreadSynchronize();

    crefTKernel<<<blocksPerGrid, threadsPerBlock>>>(
            devRefTData,
            devCrefTOffsets,
            devTgtData,
            devBHData,
            devCrefTData,
            devTgtTmpData,
            devCount0Data,
            (int)(tgtData.size()));

    cudaThreadSynchronize();
    checkCUDAError(cudaGetLastError())

    threadsPerBlock = pbwtThreadsPerTgt;
    // NOTE: Kwords refers to sizeof(data_type) words. We need the same space as for the crefT data times 2 for one (fwd or bck) compressed PBWT. x4 is because we count uint32_t words now
    uint32_t* devPBWTData;
    checkCUDAError(cudaMalloc(&devPBWTData, 4*crefTOffsets.back()*Kwords*sizeof(uint32_t)))

    cudaThreadSynchronize();

//    // DEBUG
//    cout << "GPU PBWT kernel (fwd): " << endl;
//    cout << " Threads per block: " << threadsPerBlock << endl;
//    cout << " Blocks per grid:   " << blocksPerGrid << endl;

    // forward PBWT
    pbwtPartedKernel<<<blocksPerGrid, threadsPerBlock>>>(
            devCrefTOffsets,
            devCrefTData,
            devPBWTData,
            (uint32_t*) devTgtTmpData,
            (int)(tgtData.size()),
            false); // forward

    cudaThreadSynchronize();
    checkCUDAError(cudaGetLastError())

    // copy result data to host
    for (size_t t = 0; t < tgtData.size(); t++) {
//        cout << t << "/" << tgtData.size() << ": " << hex << (size_t) crefT[t] << " <- " << (size_t) (devCrefTData+crefTOffsets[t]) << ": " << dec << crefTOffsets[t] << " - " << crefTOffsets[t+1] << endl;
        // NOTE: Kwords refers to sizeof(data_type), but here, we copy sizeof(uint32_t). Due to forward PBWT + offsets, we need to multiply with 4
        checkCUDAError(cudaMemcpy(cPBWT[t], devPBWTData + crefTOffsets[t]*Kwords*4, (crefTOffsets[t+1]-crefTOffsets[t])*Kwords*4*sizeof(uint32_t), cudaMemcpyDeviceToHost))
        checkCUDAError(cudaMemcpy(count0s[t]->data(), devCount0Data+crefTOffsets[t], (crefTOffsets[t+1]-crefTOffsets[t])*sizeof(int), cudaMemcpyDeviceToHost))
    }

    cudaThreadSynchronize();

//    // DEBUG
//    cout << "GPU PBWT kernel (bck). " << endl;

    // backward PBWT
    pbwtPartedKernel<<<blocksPerGrid, threadsPerBlock>>>(
            devCrefTOffsets,
            devCrefTData,
            devPBWTData,
            (uint32_t*) devTgtTmpData,
            (int)(tgtData.size()),
            true); // reverse

    cudaThreadSynchronize();
    checkCUDAError(cudaGetLastError())

    // copy result data to host
    for (size_t t = 0; t < tgtData.size(); t++) {
//        cout << t << "/" << tgtData.size() << ": " << hex << (size_t) crefT[t] << " <- " << (size_t) (devCrefTData+crefTOffsets[t]) << ": " << dec << crefTOffsets[t] << " - " << crefTOffsets[t+1] << endl;
        // NOTE: Kwords refers to sizeof(data_type), but here, we copy sizeof(uint32_t). Due to backward PBWT + offsets, we need to multiply with 4 and need an offset at the destination to jump over forward PBWT
        checkCUDAError(cudaMemcpy(cPBWT[t] + (crefTOffsets[t+1]-crefTOffsets[t])*Kwords*4, devPBWTData + crefTOffsets[t]*Kwords*4, (crefTOffsets[t+1]-crefTOffsets[t])*Kwords*4*sizeof(uint32_t), cudaMemcpyDeviceToHost))
    }

    // free memory for return data
    cudaFree(devCrefTData);
    cudaFree(devPBWTData);
    cudaFree(devCount0Data);
}

void GPUEngine::free() {
    cudaFree(devRefTData);
    cudaFree(devCrefTOffsets);
    cudaFree(devTgtData);
    cudaFree(devBHData);
    cudaFree(devTgtTmpData);
}

#endif // USE_CUDA_GPU
