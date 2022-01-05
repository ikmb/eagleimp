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

#ifndef GPUHANDLER_H_
#define GPUHANDLER_H_

#ifdef USE_CUDA_GPU

#include <atomic>
#include <memory>
#include <chrono>
#include <iostream>
#include <utility>
#include <time.h>

#include "hybridsys/Hybridsys.h"
#include "hybridsys/ThreadUtils.h"
#include "hybridsys/Buffer.h"

#include "Target.h"
#include "FPGAConfigurationEagleImp.h"
#include "FPGAHandler.h"
#include "GPUEngine.h"
#include "ThreadPool.h"
#include "VCFData.h"
#include "PBWTPhaser.h"

using namespace std;

class GPUHandler {

public:

    GPUHandler(hybridsys::Hybridsys &hysys, const VCFData &vcfdata);

    // initializes all GPU engines and returns the block size (i.e. the number of targets per block)
    size_t initialize(size_t K, size_t Nrefhapscorr, size_t NrefhapsAD);

    void preProcessGPUOnly(
                tbb::concurrent_bounded_queue<PBWTPhaser::targetQueue_type> &inqueue,
                tbb::concurrent_bounded_queue<PBWTPhaser::cPBWTQueue_type> &outqueue,
                atomic_bool& termination_request,
                int threadIndex
                );

    void processCondensedRef(
                tbb::concurrent_bounded_queue<PBWTPhaser::cPBWTQueue_type> &inqueue,
                tbb::concurrent_bounded_queue<PBWTPhaser::confidence_type> &outqueue,
                atomic_bool& termination_request,
                int threadIndex
                );


private:
    hybridsys::Hybridsys &hysys;
    std::vector<GPUEngine> engines;
    const VCFData &vcfdata;

    size_t M;
    size_t K = 0;
    size_t Kwords = 0; // number of crefwords_type words to fit K
    size_t Nrefhapscorr = 0;
    size_t NrefhapsAD= 0;
//    int reserveTargets;
//    int maxTargets;
//    int maxM;

    size_t totalTargets;
    atomic<size_t> targetsRemainingPreProcess;
    atomic<size_t> targetsRemainingProcessCref;
    size_t totalBlocks = 0;
    size_t gpuBlockSize = 0; // the number of targets to be analyzed with one kernel call
};

#endif // USE_CUDA_GPU

#endif /* GPUHANDLER_H_ */
