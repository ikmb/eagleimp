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

#ifndef FPGAHANDLER_H_
#define FPGAHANDLER_H_

#ifdef USE_AD_FPGA

#include <atomic>
#include <memory>
#include <chrono>
#include <iostream>
#include <utility>
#include <time.h>
#include <omp.h>

#include "hybridsys/Hybridsys.h"
#include "hybridsys/ThreadUtils.h"
#include "hybridsys/Buffer.h"

#include "Target.h"
#include "FPGAConfigurationEagleImp.h"
#include "ThreadPool.h"
#include "VCFData.h"
#include "PBWTPhaser.h"

using namespace std;

class FPGAHandler
{
public:


//    typedef vector<BooleanVector> tgtsplitsites_type;
//    typedef vector<uint32_t> tgtnsplitsites_type;

    FPGAHandler(
            hybridsys::Hybridsys &hysys,
            chrono::milliseconds timeout,
            hybridsys::BufferFactory<hybridsys::FPGABuffer> &fpgafactory_wr,
            hybridsys::BufferFactory<hybridsys::FPGABuffer> &fpgafactory_rd,
            hybridsys::BufferFactory<hybridsys::CUDABuffer> &gpufactory,
            const VCFData &vcfdata,
            size_t maxsplitsites,
            unsigned numthreads,
            bool fpgaOnly,
            bool debug = false);

    ~FPGAHandler() {
//        free(dbgbufs[0]);
//        free(dbgbufs[1]);
//        free(dbgbufs[2]);
    }

    // initializes next iteration
    void initIteration(size_t K, size_t iter);

//    // analyses targets in the block and provides the split sites as a BooleanVector and the count for each target
//    void findSplitSites(size_t firsttarget, vector<BooleanVector> &splitsites, vector<uint32_t> &nsplitsites);

    size_t getBlockSize() const { return blocksize; }
    size_t getTotalBlocks() const { return totalBlocks; }

    // initializes transmission buffers with reference data, re-used for every target block
    void prepareReferenceBuffers();

    // initializes constants and target data from first_target to first_target+fpgaconf.getNumPipelines() (only one target if FPGAs are not used)
    void initTargets(const vector<Target*>&, bool lastblock, int threadIndex);

    // processes a target block by calling initTargets and sending the references accordingly
    void processBlock(const vector<Target*> &target_block, bool lastblock, int threadIndex);

    // thread function that prepares the targets to be processed on the FPGA
    void preparePhasing(
            tbb::concurrent_bounded_queue<Target*> &inqueue,
            tbb::concurrent_bounded_queue<Target*> &outqueue,
            atomic_bool& termination_request,
            int threadIndex
            );

    // thread function that initializes the targets on the FPGA and streams the references
    void provideReferences(
            tbb::concurrent_bounded_queue<Target*> &inqueue,
            tbb::concurrent_bounded_queue<PBWTPhaser::targetQueue_type> &outqueue,
            atomic_bool& termination_request,
            int threadIndex
            );

//    void readFPGA(
//            atomic_bool& termination_request,
//            int threadIndex
//            );

    // processes the returned FPGA buffers with condensed PBWT data, providing it to phaser
    void preProcessFPGAOnly(
            tbb::concurrent_bounded_queue<PBWTPhaser::targetQueue_type> &inqueue,
            tbb::concurrent_bounded_queue<PBWTPhaser::cPBWTQueue_type> &outqueue,
            atomic_bool& termination_request,
            int threadIndex
            );

    // thread function that processes the PBWT from the FPGA and phases the targets, returns the phase confidences in the outqueue
    void processPBWT(
            tbb::concurrent_bounded_queue<PBWTPhaser::cPBWTQueue_type> &inqueue,
            tbb::concurrent_bounded_queue<PBWTPhaser::confidence_type> &outqueue,
            atomic_bool& termination_request,
            int threadIndex
            );

    bool isFinished() {
        return targetsRemaining_out == 0;
    }

    // shows the progress of the FPGA buffers between 0 (not started) and 1 (finished)
    double progress() {
        return 1.0 - ((double)targetsRemaining_out / totalTargets);
    }

    // how many uint32_t words are required for general initialization (including sync word, excluding target specific data)
    // here: 256 bit start word + 256 bit constants
    static constexpr int num_constants = 16;

private:
    hybridsys::Hybridsys &hysys;
    FPGAConfigurationEagleImp fpgaconf;
    chrono::milliseconds timeout;
    hybridsys::BufferFactory<hybridsys::FPGABuffer> &bufferFactoryFPGA_write;
    hybridsys::BufferFactory<hybridsys::FPGABuffer> &bufferFactoryFPGA_read;
    hybridsys::BufferFactory<hybridsys::CUDABuffer> &bufferFactoryGPU; // only required if used in combination with GPU
    const VCFData &vcfdata;
    vector<hybridsys::FPGABuffer> refbuffers;
    size_t maxpbwtsites;
    size_t K = 0;
    size_t iter = 0;

//    tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>> fpgaReadQueue;

    size_t totalTargets;
    atomic<size_t> targetsRemaining_in;
    atomic<size_t> targetsRemaining_out;
    size_t totalBlocks;
    size_t blocksize;
    atomic<size_t> fpga_blocks_sent;
    atomic<bool> last_block_flag;

    size_t constrequired_bufsize;
    size_t tgtsize;

    unsigned numthreads;

    bool fpgaOnly;

    bool debug;

//    array<char*,3> dbgbufs;
};

#endif // USE_AD_FPGA

#endif /* FPGAHANDLER_H_ */
