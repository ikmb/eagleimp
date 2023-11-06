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

#include <iostream>
#include <iomanip>
#include <omp.h>

#include "hybridsys/Hybridsys.h"

#include "FPGAHandler.h"
#include "GPUHandler.h"
#include "Target.h"
#include "ThreadPool.h"
#include "HaplotypePath.h"
#include "utils.h"
#include "StatusFile.h"
#include "Stopwatch.h"

#include "PBWTPhaser.h"

using namespace std;
using namespace placeholders;
using namespace hybridsys;

/* static */ atomic<bool> PBWTPhaser::terminate(false);

PBWTPhaser::PBWTPhaser(Hybridsys &hysys_, VCFData &vcfdata_, unsigned numthreads_, size_t Karg_, uint32_t iters_,
        fp_type expectIBDcM_, fp_type hist_, fp_type pErr_, fp_type pLimit_, bool impMissing_, bool doPrePhasing_, bool noRevPhasing_, bool skipPhasing_, bool debug_)
    :
      numthreads(numthreads_),
      vcfdata(vcfdata_),
//      targetsFull(vcfdata_.getTargets()),
      nTarget(vcfdata_.getNTarget()),
      Karg(Karg_),
      iters(iters_),
      expectIBDcM(expectIBDcM_),
      hist(hist_),
      pErr(pErr_),
      pLimit(pLimit_),
      hysys(hysys_),
      impMissing(impMissing_),
      doPrePhasing(doPrePhasing_),
      doRevPhasing(!noRevPhasing_),
      skipPhasing(skipPhasing_),
      targetIDs(vcfdata.getTargetIDs()),
      debug(debug_)
{}

void PBWTPhaser::phase(vector<BooleanVector> &phasedTargets, vector<vector<float>> &phasedDosages, vector<fp_type> &totconfidences, vector<size_t> &ncalls, int chunk) {

    atomic<int> missings(0);
    atomic<int> monhets(0);

    // do the phasing

    if (usefpga) {

        phaseFPGA(phasedTargets, phasedDosages, chunk, totconfidences, ncalls);

    } else if (usegpu) {

        phaseGPU(phasedTargets, phasedDosages, chunk, totconfidences, ncalls);

    } else { // no FPGA/GPU acceleration

        phaseCPU(phasedTargets, phasedDosages, chunk, totconfidences, ncalls);

    } // END no FPGA/GPU

}

void PBWTPhaser::phaseFPGA(vector<BooleanVector> &phasedTargets __attribute__((unused)), vector<vector<float>> &phasedDosages __attribute__((unused)), int chunk __attribute__((unused)),
        vector<fp_type> &totconfidences __attribute__((unused)), vector<size_t> &ncalls __attribute__((unused))) {

#ifdef USE_AD_FPGA

    cout << "Using FPGA support and " <<
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT //|| defined STOPWATCH
             1
#else
            numthreads
#endif
            << " threads for phasing." << endl;

    // init FPGA buffers
    BufferFactory<FPGABuffer> fpgaBufferFactory(fpga_buffersize);
    for(unsigned j = 0; j < fpga_num_buffers; ++j) {
        fpgaBufferFactory.preallocateBuffer();
    }

    // init GPU buffers
//#ifdef USE_CUDA_GPU
    BufferFactory<CUDABuffer> gpuBufferFactory(gpu_buffersize);
//#else
//    BufferFactory<FPGABuffer> gpuBufferFactory(gpu_buffersize);
//#endif
    if (usegpu) {
        for(unsigned j = 0; j < gpu_num_buffers; ++j) {
            gpuBufferFactory.preallocateBuffer();
        }
    }

    // set up transmission queues
    tbb::concurrent_bounded_queue<targetQueue_type> provideQueue;
    tbb::concurrent_bounded_queue<targetQueue_type> processQueue;
//    tbb::concurrent_bounded_queue<condensedRefQueue_type> crefFPGAQueue;
    tbb::concurrent_bounded_queue<cPBWTQueue_type> cPBWTQueue;
    tbb::concurrent_bounded_queue<confidence_type> confidenceQueue; // contains only target ID and corresponding phasing confidence (incl. number of call sites) for result printing
    cPBWTQueue.set_capacity(2*numthreads); // if this is larger and the FPGA is faster than the phasing process, we can quickly get out of memory here for large datasets!
    confidenceQueue.set_capacity(nTarget);

    stringstream ss;
    ss << "Phasing (Chunk " << chunk+1 << "/" << vcfdata.getNChunks() << ")";
    string statusstring = ss.str();
    StatusFile::updateStatus(0, statusstring);

    FPGAHandler fpgahandler(hysys, fpga_timeout, fpgaBufferFactory, gpuBufferFactory, vcfdata, maxpbwtsites, numthreads, !usegpu, debug);

    for (uint32_t iter = 1; iter <= iters; iter++) {

        // the actual K may not be bigger than the number of provided haplotypes, but has to be adjusted according to the current phasing iteration
        size_t localNrefhapsAD;
        size_t localNrefhapsCorr;
        size_t localK;
        // TODO triple(!!!) code! -> make function!
        if (iter > 1) {
            localNrefhapsAD = vcfdata.getNReferenceHapsMax(); // NOTE: already corrected when DEBUG_TARGET* set
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
            localNrefhapsCorr = localNrefhapsAD - vcfdata.getNChrXYHaploidsRef() -2; // NOTE: from the available haps we need to remove the one diploid target (since it does not use it's own phasing result as reference in further iterations)
            for (int dbgtgt = DEBUG_TARGET_START; dbgtgt <= DEBUG_TARGET_STOP; dbgtgt++) {
                if (vcfdata.getChrXYHaploidsTgt()[dbgtgt])
                    localNrefhapsCorr--; // reduce one hap for each haploid target in debug range
            }
#else
            localNrefhapsCorr = localNrefhapsAD - vcfdata.getNHaploidsRef() - vcfdata.getNHaploidsTgt() -2; // NOTE: from the available haps we need to remove the one diploid target (since it does not use it's own phasing result as reference in further iterations)
#endif

            // add transposition of phased targets to transposed reference for upcoming iteration
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
            size_t targetcapacity = phasedTargets[DEBUG_TARGET_START].getCapacity();
            size_t targetsize = phasedTargets[DEBUG_TARGET_START].size();
#else
            size_t targetcapacity = phasedTargets[0].getCapacity();
            size_t targetsize = phasedTargets[0].size();
#endif
            size_t NrefhapsAD = vcfdata.getNReferenceHaps();
            size_t NrefhapsADwords = divideRounded(NrefhapsAD, 8*sizeof(BooleanVector::data_type));
            size_t hapoff = NrefhapsAD % (8*sizeof(BooleanVector::data_type));
            size_t mask = hapoff == 0 ? 0ull : ((1ull << hapoff)-1);
            size_t clearsize = divideRounded(2*nTarget - (hapoff ? (8*sizeof(BooleanVector::data_type) - hapoff) : 0), (size_t) 8);
            for (size_t m = 0; m < targetsize; m++) {
                // if this is the third or further iteration, we need to clear the space first
                if (iter > 2) {
                    vcfdata.getReferenceT()[m].getData()[NrefhapsADwords-1] &= mask;
                    memset(vcfdata.getReferenceT()[m].getData()+NrefhapsADwords, 0, clearsize);
                }
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                for (size_t nt = DEBUG_TARGET_START; nt <= DEBUG_TARGET_STOP; nt++) {
                    vcfdata.getReferenceT()[m].setPairWithPreInit(NrefhapsAD + 2*nt - 2*DEBUG_TARGET_START, phasedTargets[2*nt][m], phasedTargets[2*nt+1][m]);
                }
                vcfdata.getReferenceT()[m].setSize(NrefhapsAD + 2*nTarget - 2*DEBUG_TARGET_START);
#else
                for (size_t nt = 0; nt < nTarget; nt++) {
                    vcfdata.getReferenceT()[m].setPairWithPreInit(NrefhapsAD + 2*nt, phasedTargets[2*nt][m], phasedTargets[2*nt+1][m]);
                }
                vcfdata.getReferenceT()[m].setSize(NrefhapsAD + 2*nTarget);
#endif
            }

            // append phased targets to reference for upcoming iteration
            vector<BooleanVector>::iterator refit = vcfdata.getReference().begin();
            refit += vcfdata.getNReferenceHaps(); // jump to beginning of target space
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
            for (size_t nthap = 2*DEBUG_TARGET_START; nthap <= 2*DEBUG_TARGET_STOP+1; nthap++, refit++)
#else
            for (size_t nthap = 0; nthap < 2*nTarget; nthap++, refit++)
#endif
            {
                memcpy(refit->getData(), phasedTargets[nthap].getData(), targetcapacity);
                refit->setSize(targetsize);
                phasedTargets[nthap].clear();
            }

        } else { // iter == 0
            localNrefhapsAD = vcfdata.getNReferenceHaps();
            localNrefhapsCorr = localNrefhapsAD - vcfdata.getNHaploidsRef();
        }
        localK = min(localNrefhapsCorr, Karg);
        // FPGA supports only even K
        if (localK % 2)
            localK--;
        // acquire FPGA lock
        StatusFile::updateStatus(0, statusstring + ": Waiting");
        fpgahandler.initIteration(localK, iter);
        // if we received a user termination request during waiting for the FPGA lock, we can stop here
        if (terminate) {
            for(FPGA &f : hysys.getFPGAs())
                f.unlock();
            cerr << "User abort." << endl;
            exit(EXIT_FAILURE);
        }
        StatusFile::updateStatus(0, statusstring);

        fpgahandler.prepareReferenceBuffers();

        cout << "Phasing iteration " << iter << "/" << iters << " (K=" << localK << "): 0%" << flush;
        int pgb = 0; // for progress bar

        // processor (currently only 1 thread) that initializes the targets to the FPGA and streams the references
        ThreadPool<targetQueue_type, targetQueue_type> fpgaProvider(
                1,
                provideQueue,
                processQueue,
                bind(&FPGAHandler::provideReferences, &fpgahandler, _1, _2, _3, _4));

        // processor (currently only 1 thread) that collects the condensed references created by the FPGA and prepares
        // targets that can be phased multi-threaded by the fpgaPhasingProcessor
        ThreadPool<targetQueue_type, cPBWTQueue_type> fpgaOnlyPreProcessor(
                1,
                processQueue,
                cPBWTQueue,
                bind(&FPGAHandler::preProcessFPGAOnly, &fpgahandler, _1, _2, _3, _4));

        // processor for phasing the targets where the FPGA prepared the condensed references for, with PBWT generation on the host
        ThreadPool<cPBWTQueue_type, confidence_type> fpgaPhasingProcessor(
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT || defined STOPWATCH
                1,
#else
                numthreads, // use user assigned number of threads
#endif
                cPBWTQueue,
                confidenceQueue,
                bind(&FPGAHandler::processPBWT, &fpgahandler, _1, _2, _3, _4));

        // start phasing thread using the condensed references from the FPGA
        fpgaPhasingProcessor.run();
        // start fetching FPGA data and preprocessing for FPGA only
        fpgaOnlyPreProcessor.run();
        // start providing target/reference data to FPGA
        fpgaProvider.run();

        // fill inqueue of provide thread
        size_t curr_filled = 0;
        size_t blocksize = fpgahandler.getBlockSize();
        vector<Target*> *target_block = new vector<Target*>;
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
        for (size_t nt = DEBUG_TARGET_START; nt <= DEBUG_TARGET_STOP; nt++)
#else
        for (size_t nt = 0; nt < nTarget && !terminate; nt++)
#endif
        {
            if (curr_filled == blocksize) { // need to push the current data to the queue
                if (debug)
                    cout << "Pushing block. (curr target: " << nt << ")" << endl;

                provideQueue.push(shared_ptr<vector<Target*>>(target_block));

                target_block = new vector<Target*>; // assign new block
                curr_filled = 0;
            }
            Target* tgtptr = new Target(phasedTargets[2*nt], phasedTargets[2*nt+1], phasedDosages[2*nt], phasedDosages[2*nt+1],
                    nt, targetIDs[nt], vcfdata, localNrefhapsAD, localNrefhapsCorr, localK, expectIBDcM, hist, pErr, pLimit, doPrePhasing,
                    doRevPhasing, iter == iters, impMissing, skipPhasing, usefpga);
            target_block->push_back(tgtptr);
            curr_filled++;
        }

        // push the last data (probably containing less than blocksize targets, but at least one!)
        if (debug)
            cout << "Pushing last block. (size: " << target_block->size() << ")" << endl;
        // sets the termination request and ensures that at least one block is in the queue to release a pop operation
        // and to mark the last block transmitted to the FPGA as "last_block".
        if (terminate) {
            fpgaProvider.setTerminationRequest(false); // do not abort inqueue! the next pop operation has to be successful to mark the block sent to the FPGA as "last_block"
            fpgaOnlyPreProcessor.setTerminationRequest(false); // do not abort inqueue as the targets in the queue were all going to be processed by the FPGA and need to be fetched!
            fpgaPhasingProcessor.setTerminationRequest(true); // abort inqueue as the preliminary process will not push to its outqueue after receiving the termination request
        }
        provideQueue.push(shared_ptr<vector<Target*>>(target_block));

        // just fetch all confidences
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
        for (size_t nt = DEBUG_TARGET_START; nt <= DEBUG_TARGET_STOP; nt++)
#else
        for (size_t nt = 0; nt < nTarget; nt++)
#endif
        { // we might not know the order in which the results arrive, but we know that there must be exactly nTarget results
            if (nt % 64 == 0) { // update status file every 64 targets
//                float progress = ((nt/(float)nTarget)/(float)iters)/(float)vcfdata.getNChunks();
//                progress += ((iter-1)/(float)iters)/(float)vcfdata.getNChunks();
//                progress += chunk/(float)vcfdata.getNChunks();
                StatusFile::updateStatus(nt/(float)nTarget);
                if (pgb==3) { // printing % after three dots
                    cout << (100*nt)/nTarget << "%" << flush;
                    pgb = 0;
                } else {
                    cout << "." << flush;
                    pgb++;
                }
            }
            // terminate the processor threads if there's a user termination request
            if (terminate) { // termination request received -> terminate processor threads and leave
                fpgaProvider.setTerminationRequest(false); // do not abort inqueue! the next pop operation has to be successful to mark the block sent to the FPGA as "last_block"
                fpgaOnlyPreProcessor.setTerminationRequest(false); // do not abort inqueue as the targets in the queue were all going to be processed by the FPGA and need to be fetched!
                fpgaPhasingProcessor.setTerminationRequest(true); // abort inqueue as the preliminary process will not push to its outqueue after receiving the termination request
                break;
            }
            confidence_type tconf;
            confidenceQueue.pop(tconf);
            totconfidences[tconf.id] += tconf.totalconf;
            ncalls[tconf.id] += tconf.ncalls;
        }

        if (terminate) { // we could terminate now
            fpgaProvider.waitForTermination();
            fpgaOnlyPreProcessor.waitForTermination();
            fpgaPhasingProcessor.waitForTermination();
            cerr << "\nUser abort." << endl;
            exit(EXIT_FAILURE);
        }

        if (pgb==0) // just printed "xx%"
            cout << ".";
        cout << "100%" << endl;

        StatusFile::nextStep(); // each iteration is one step!
    }

    // if we manipulated the reference by adding the phased targets, we need to rewind this change at least for the transposed reference
    if (iters > 1) {
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
        size_t targetsize = phasedTargets[DEBUG_TARGET_START].size();
#else
        size_t targetsize = phasedTargets[0].size();
#endif
        size_t NrefhapsAD = vcfdata.getNReferenceHaps();
        size_t NrefhapsADwords = divideRounded(NrefhapsAD, 8*sizeof(BooleanVector::data_type));
        size_t hapoff = NrefhapsAD % (8*sizeof(BooleanVector::data_type));
        size_t mask = hapoff == 0 ? 0ull : ((1ull << hapoff)-1);
        size_t clearsize = divideRounded(2*nTarget - (hapoff ? (8*sizeof(BooleanVector::data_type) - hapoff) : 0), (size_t) 8);
        for (size_t m = 0; m < targetsize; m++) {
            // clear the space and set the length
            vcfdata.getReferenceT()[m].getData()[NrefhapsADwords-1] &= mask;
            memset(vcfdata.getReferenceT()[m].getData()+NrefhapsADwords, 0, clearsize);
            vcfdata.getReferenceT()[m].setSize(NrefhapsAD);
        }
    }
#endif // USE_AD_FPGA
}

void PBWTPhaser::phaseGPU(vector<BooleanVector> &phasedTargets __attribute__((unused)), vector<vector<float>> &phasedDosages __attribute__((unused)), int chunk __attribute__((unused)),
        vector<fp_type> &totconfidences __attribute__((unused)), vector<size_t> &ncalls __attribute__((unused))) {

#ifdef USE_CUDA_GPU

    {
        stringstream ss;
        ss << "Phasing (Chunk " << chunk+1 << "/" << vcfdata.getNChunks() << ")";
        StatusFile::updateStatus(chunk/(float)vcfdata.getNChunks(), ss.str());
    }

//    BufferFactory<CUDABuffer> gpuBufferFactory(gpu_buffersize);
//    for(unsigned j = 0; j < gpu_num_buffers; ++j) {
//        gpuBufferFactory.preallocateBuffer();
//    }
//    cout << "Preloaded GPU buffers." << endl;

    cout << "Using GPU support and " <<
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT //|| defined STOPWATCH
             1
#else
            numthreads
#endif
            << " threads for phasing." << endl;

    // set up transmission queues
    tbb::concurrent_bounded_queue<targetQueue_type> provideQueue;
    tbb::concurrent_bounded_queue<cPBWTQueue_type> cpbwtGPUQueue;
    tbb::concurrent_bounded_queue<confidence_type> confidenceQueue; // contains only target ID and corresponding phasing confidence (incl. number of call sites) for result printing
    confidenceQueue.set_capacity(nTarget);

    GPUHandler gpuhandler(hysys, vcfdata);

    for (uint32_t iter = 1; iter <= iters; iter++) {
        Stopwatch swpi(string("iter ").append(to_string(iter)).c_str());

        // the actual K may not be bigger than the number of provided haplotypes, but has to be adjusted according to the current phasing iteration
        size_t localNrefhapsAD;
        size_t localNrefhapsCorr;
        size_t localK;
        // TODO triple(!!!) code! -> make function!
        if (iter > 1) {
            localNrefhapsAD = vcfdata.getNReferenceHapsMax(); // NOTE: already corrected when DEBUG_TARGET* set
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
            localNrefhapsCorr = localNrefhapsAD - vcfdata.getNChrXYHaploidsRef() -2; // NOTE: from the available haps we need to remove the one diploid target (since it does not use it's own phasing result as reference in further iterations)
            for (int dbgtgt = DEBUG_TARGET_START; dbgtgt <= DEBUG_TARGET_STOP; dbgtgt++) {
                if (vcfdata.getChrXYHaploidsTgt()[dbgtgt])
                    localNrefhapsCorr--; // reduce one hap for each haploid target in debug range
            }
#else
            localNrefhapsCorr = localNrefhapsAD - vcfdata.getNChrXYHaploidsRef() - vcfdata.getNChrXYHaploidsTgt() -2; // NOTE: from the available haps we need to remove the one diploid target (since it does not use it's own phasing result as reference in further iterations)
#endif

            // add transposition of phased targets to transposed reference for upcoming iteration
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
            size_t targetcapacity = phasedTargets[DEBUG_TARGET_START].getCapacity();
            size_t targetsize = phasedTargets[DEBUG_TARGET_START].size();
#else
            size_t targetcapacity = phasedTargets[0].getCapacity();
            size_t targetsize = phasedTargets[0].size();
#endif
            size_t NrefhapsAD = vcfdata.getNReferenceHaps();
            size_t NrefhapsADwords = divideRounded(NrefhapsAD, 8*sizeof(BooleanVector::data_type));
            size_t hapoff = NrefhapsAD % (8*sizeof(BooleanVector::data_type));
            size_t mask = hapoff == 0 ? 0ull : ((1ull << hapoff)-1);
            size_t clearsize = divideRounded(2*nTarget - (hapoff ? (8*sizeof(BooleanVector::data_type) - hapoff) : 0), (size_t) 8);
            for (size_t m = 0; m < targetsize; m++) {
                // if this is the third or further iteration, we need to clear the space first
                if (iter > 2) {
                    vcfdata.getReferenceT()[m].getData()[NrefhapsADwords-1] &= mask;
                    memset(vcfdata.getReferenceT()[m].getData()+NrefhapsADwords, 0, clearsize);
                }
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                for (size_t nt = DEBUG_TARGET_START; nt <= DEBUG_TARGET_STOP; nt++) {
                    vcfdata.getReferenceT()[m].setPairWithPreInit(NrefhapsAD + 2*nt - 2*DEBUG_TARGET_START, phasedTargets[2*nt][m], phasedTargets[2*nt+1][m]);
                }
                vcfdata.getReferenceT()[m].setSize(NrefhapsAD + 2*nTarget - 2*DEBUG_TARGET_START);
#else
                for (size_t nt = 0; nt < nTarget; nt++) {
                    vcfdata.getReferenceT()[m].setPairWithPreInit(NrefhapsAD + 2*nt, phasedTargets[2*nt][m], phasedTargets[2*nt+1][m]);
                }
                vcfdata.getReferenceT()[m].setSize(NrefhapsAD + 2*nTarget);
#endif
            }

            // append phased targets to reference for upcoming iteration
            vector<BooleanVector>::iterator refit = vcfdata.getReference().begin();
            refit += vcfdata.getNReferenceHaps(); // jump to beginning of target space
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
            for (size_t nthap = 2*DEBUG_TARGET_START; nthap <= 2*DEBUG_TARGET_STOP+1; nthap++, refit++)
#else
            for (size_t nthap = 0; nthap < 2*nTarget; nthap++, refit++)
#endif
            {
                memcpy(refit->getData(), phasedTargets[nthap].getData(), targetcapacity);
                refit->setSize(targetsize);
                phasedTargets[nthap].clear();
            }

        } else { // iter == 0
            localNrefhapsAD = vcfdata.getNReferenceHaps();
            localNrefhapsCorr = localNrefhapsAD - vcfdata.getNChrXYHaploidsRef();
        }
        localK = min(localNrefhapsCorr, Karg);

        size_t blocksize = gpuhandler.initialize(localK, localNrefhapsCorr, localNrefhapsAD);

        cout << "Phasing iteration " << iter << "/" << iters << " (K=" << localK << "): 0%" << flush;
        int pgb = 0; // for progress bar

        // processor that provides the targets to the GPU, starts GPU kernels for creating the condensed references on the GPU, and copies the sequences back to the host
        ThreadPool<targetQueue_type, cPBWTQueue_type> gpuOnlyPreProcessor(
                hysys.getGPUs().size(),
                provideQueue,
                cpbwtGPUQueue,
                bind(&GPUHandler::preProcessGPUOnly, &gpuhandler, _1, _2, _3, _4));

        // processor for phasing the targets where the GPU prepared the condensed references for, with PBWT generation on the host
        ThreadPool<cPBWTQueue_type, confidence_type> gpuPhasingProcessor(
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT || defined STOPWATCH
                1,
#else
                numthreads, // use user assigned number of threads
#endif
                cpbwtGPUQueue,
                confidenceQueue,
                bind(&GPUHandler::processCondensedRef, &gpuhandler, _1, _2, _3, _4));

        gpuPhasingProcessor.run();
        gpuOnlyPreProcessor.run();

        // fill inqueue of preprocessor thread
        size_t curr_filled = 0;
        vector<Target*> *target_block = new vector<Target*>;
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
        for (size_t nt = DEBUG_TARGET_START; nt <= DEBUG_TARGET_STOP; nt++)
#else
        for (size_t nt = 0; nt < nTarget; nt++)
#endif
        { // TODO implement user abort
            if (curr_filled == blocksize) { // need to push the current data to the queue
                if (debug)
                    cout << "Pushing block. (curr target: " << nt << ")" << endl;

                provideQueue.push(shared_ptr<vector<Target*>>(target_block));

                target_block = new vector<Target*>; // assign new block
                curr_filled = 0;
            }
            Target* tgtptr = new Target(phasedTargets[2*nt], phasedTargets[2*nt+1], phasedDosages[2*nt], phasedDosages[2*nt+1],
                    nt, targetIDs[nt], vcfdata, localNrefhapsAD, localNrefhapsCorr, localK, expectIBDcM, hist, pErr, pLimit, doPrePhasing,
                    doRevPhasing, iter == iters, impMissing, skipPhasing, usegpu);
            target_block->push_back(tgtptr);
            curr_filled++;
        }
        // push the last data (probably containing less than blocksize targets)
        if (debug)
            cout << "Pushing last block. (size: " << target_block->size() << ")" << endl;

        blocksize = target_block->size(); // only required for the workaround below
        provideQueue.push(shared_ptr<vector<Target*>>(target_block));

        // just fetch all confidences
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
        for (size_t nt = DEBUG_TARGET_START; nt <= DEBUG_TARGET_STOP; nt++)
#else
        for (size_t nt = 0; nt < nTarget; nt++)
#endif
        { // we might not know the order in which the results arrive, but we know that there must be exactly nTarget results
            if (nt % 64 == 0) { // update status file every 64 targets
                float progress = ((nt/(float)nTarget)/(float)iters)/(float)vcfdata.getNChunks();
                progress += ((iter-1)/(float)iters)/(float)vcfdata.getNChunks();
                progress += chunk/(float)vcfdata.getNChunks();
                StatusFile::updateStatus(progress);
                if (pgb==3) { // printing % after three dots
                    cout << (100*nt)/nTarget << "%" << flush;
                    pgb = 0;
                } else {
                    cout << "." << flush;
                    pgb++;
                }
            }
            confidence_type tconf;
            confidenceQueue.pop(tconf);
            totconfidences[tconf.id] += tconf.totalconf;
            ncalls[tconf.id] += tconf.ncalls;
        }
        if (pgb==0) // just printed "xx%"
            cout << ".";
        cout << "100%" << endl;
        swpi.stop();

        StatusFile::nextStep(); // each iteration is one step!
    }

    // if we manipulated the reference by adding the phased targets, we need to rewind this change at least for the transposed reference
    if (iters > 1) {
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
        size_t targetsize = phasedTargets[DEBUG_TARGET_START].size();
#else
        size_t targetsize = phasedTargets[0].size();
#endif
        size_t NrefhapsAD = vcfdata.getNReferenceHaps();
        size_t NrefhapsADwords = divideRounded(NrefhapsAD, 8*sizeof(BooleanVector::data_type));
        size_t hapoff = NrefhapsAD % (8*sizeof(BooleanVector::data_type));
        size_t mask = hapoff == 0 ? 0ull : ((1ull << hapoff)-1);
        size_t clearsize = divideRounded(2*nTarget - (hapoff ? (8*sizeof(BooleanVector::data_type) - hapoff) : 0), (size_t) 8);
        for (size_t m = 0; m < targetsize; m++) {
            // clear the space and set the length
            vcfdata.getReferenceT()[m].getData()[NrefhapsADwords-1] &= mask;
            memset(vcfdata.getReferenceT()[m].getData()+NrefhapsADwords, 0, clearsize);
            vcfdata.getReferenceT()[m].setSize(NrefhapsAD);
        }
    }
#endif // USE_CUDA_GPU
}

void PBWTPhaser::phaseCPU(vector<BooleanVector> &phasedTargets, vector<vector<float>> &phasedDosages, int chunk,
        vector<fp_type> &totconfidences, vector<size_t> &ncalls) {

    {
        stringstream ss;
        ss << "Phasing (Chunk " << chunk+1 << "/" << vcfdata.getNChunks() << ")";
        StatusFile::updateStatus(0, ss.str());
    }

    cout << "Using " <<
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT //|| defined STOPWATCH
                     1
#else
                    numthreads
#endif
                    << " threads for phasing." << endl;

    for (uint32_t iter = 1; iter <= iters; iter++) {
        Stopwatch swpi(string("iter ").append(to_string(iter)).c_str());

        // the actual K may not be bigger than the number of provided haplotypes, but has to be adjusted according to the current phasing iteration
        size_t localNrefhapsAD;
        size_t localNrefhapsCorr;
        size_t localK;
        // TODO triple(!!!) code! -> make function!
        if (iter > 1) {
            localNrefhapsAD = vcfdata.getNReferenceHapsMax(); // NOTE: already corrected when DEBUG_TARGET* set
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
            localNrefhapsCorr = localNrefhapsAD - vcfdata.getNChrXYHaploidsRef() -2; // NOTE: from the available haps we need to remove the one diploid target (since it does not use it's own phasing result as reference in further iterations)
            for (int dbgtgt = DEBUG_TARGET_START; dbgtgt <= DEBUG_TARGET_STOP; dbgtgt++) {
                if (vcfdata.getChrXYHaploidsTgt()[dbgtgt])
                    localNrefhapsCorr--; // reduce one hap for each haploid target in debug range
            }
#else
            localNrefhapsCorr = localNrefhapsAD - vcfdata.getNHaploidsRef() - vcfdata.getNHaploidsTgt() -2; // NOTE: from the available haps we need to remove the one diploid target (since it does not use it's own phasing result as reference in further iterations)
#endif

            // add transposition of phased targets to transposed reference for upcoming iteration
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
            size_t targetcapacity = phasedTargets[DEBUG_TARGET_START].getCapacity();
            size_t targetsize = phasedTargets[DEBUG_TARGET_START].size();
#else
            size_t targetcapacity = phasedTargets[0].getCapacity();
            size_t targetsize = phasedTargets[0].size();
#endif
            size_t NrefhapsAD = vcfdata.getNReferenceHaps();
            size_t NrefhapsADwords = divideRounded(NrefhapsAD, 8*sizeof(BooleanVector::data_type));
            size_t hapoff = NrefhapsAD % (8*sizeof(BooleanVector::data_type));
            size_t mask = hapoff == 0 ? 0ull : ((1ull << hapoff)-1);
            size_t clearsize = divideRounded(2*nTarget - (hapoff ? (8*sizeof(BooleanVector::data_type) - hapoff) : 0), (size_t) 8);
            for (size_t m = 0; m < targetsize; m++) {
                // if this is the third or further iteration, we need to clear the space first
                if (iter > 2) {
                    vcfdata.getReferenceT()[m].getData()[NrefhapsADwords-1] &= mask;
                    memset(vcfdata.getReferenceT()[m].getData()+NrefhapsADwords, 0, clearsize);
                }
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                for (size_t nt = DEBUG_TARGET_START; nt <= DEBUG_TARGET_STOP; nt++) {
                    vcfdata.getReferenceT()[m].setPairWithPreInit(NrefhapsAD + 2*nt - 2*DEBUG_TARGET_START, phasedTargets[2*nt][m], phasedTargets[2*nt+1][m]);
                }
                vcfdata.getReferenceT()[m].setSize(NrefhapsAD + 2*nTarget - 2*DEBUG_TARGET_START);
#else
                for (size_t nt = 0; nt < nTarget; nt++) {
                    vcfdata.getReferenceT()[m].setPairWithPreInit(NrefhapsAD + 2*nt, phasedTargets[2*nt][m], phasedTargets[2*nt+1][m]);
                }
                vcfdata.getReferenceT()[m].setSize(NrefhapsAD + 2*nTarget);
#endif
            }

            // append phased targets to reference for upcoming iteration
            vector<BooleanVector>::iterator refit = vcfdata.getReference().begin();
            refit += vcfdata.getNReferenceHaps(); // jump to beginning of target space
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
            for (size_t nthap = 2*DEBUG_TARGET_START; nthap <= 2*DEBUG_TARGET_STOP+1; nthap++, refit++)
#else
            for (size_t nthap = 0; nthap < 2*nTarget; nthap++, refit++)
#endif
            {
                memcpy(refit->getData(), phasedTargets[nthap].getData(), targetcapacity);
                refit->setSize(targetsize);
                phasedTargets[nthap].clear();
            }

        } else { // iter == 0
            localNrefhapsAD = vcfdata.getNReferenceHaps();
            localNrefhapsCorr = localNrefhapsAD - vcfdata.getNHaploidsRef();
        }
        localK = min(localNrefhapsCorr, Karg);

        cout << "Phasing iteration " << iter << "/" << iters << " (K=" << localK << "): 0%" << flush;
        int pgb = 0; // for progress bar

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
        for (size_t nt = DEBUG_TARGET_START; nt <= DEBUG_TARGET_STOP; nt++)
#else
#ifndef STOPWATCH
        omp_set_num_threads(numthreads);
#pragma omp parallel for schedule(dynamic,1) // chunksize 1
#endif
        for (size_t nt = 0; nt < nTarget; nt++)
#endif
        {
            if (nt % 64 == 0) { // update status file every 64 targets
//                float progress = ((nt/(float)nTarget)/(float)iters)/(float)vcfdata.getNChunks();
//                progress += ((iter-1)/(float)iters)/(float)vcfdata.getNChunks();
//                progress += chunk/(float)vcfdata.getNChunks();
                StatusFile::updateStatus(nt/(float)nTarget);
                if (pgb==3) { // printing % after three dots
                    cout << (100*nt)/nTarget << "%" << flush;
                    pgb = 0;
                } else {
                    cout << "." << flush;
                    pgb++;
                }
            }

            Target t(phasedTargets[2*nt], phasedTargets[2*nt+1], phasedDosages[2*nt], phasedDosages[2*nt+1],
                    nt, targetIDs[nt], vcfdata, localNrefhapsAD, localNrefhapsCorr, localK, expectIBDcM, hist, pErr, pLimit, doPrePhasing,
                    doRevPhasing, iter == iters, impMissing, skipPhasing, false);
            t.phase();
            totconfidences[nt] += t.getTotalConfidence();
            ncalls[nt] += t.getNCallSites();
        }
        if (pgb==0) // just printed "xx%"
            cout << ".";
        cout << "100%" << endl;
        swpi.stop();

        StatusFile::nextStep(); // each iteration is one step!
    }

    // if we manipulated the reference by adding the phased targets, we need to rewind this change at least for the transposed reference
    if (iters > 1) {
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
        size_t targetsize = phasedTargets[DEBUG_TARGET_START].size();
#else
        size_t targetsize = phasedTargets[0].size();
#endif
        size_t NrefhapsAD = vcfdata.getNReferenceHaps();
        size_t NrefhapsADwords = divideRounded(NrefhapsAD, 8*sizeof(BooleanVector::data_type));
        size_t hapoff = NrefhapsAD % (8*sizeof(BooleanVector::data_type));
        size_t mask = hapoff == 0 ? 0ull : ((1ull << hapoff)-1);
        size_t clearsize = divideRounded(2*nTarget - (hapoff ? (8*sizeof(BooleanVector::data_type) - hapoff) : 0), (size_t) 8);
        for (size_t m = 0; m < targetsize; m++) {
            // clear the space and set the length
            vcfdata.getReferenceT()[m].getData()[NrefhapsADwords-1] &= mask;
            memset(vcfdata.getReferenceT()[m].getData()+NrefhapsADwords, 0, clearsize);
            vcfdata.getReferenceT()[m].setSize(NrefhapsAD);
        }
    }

#ifdef DEBUG_PROFILE_PBWT_MAPPINGS
    {
        size_t total = 0;
        for (size_t c : HaplotypePath::dbg_mapstop)
            total += c;
        total += HaplotypePath::dbg_mapstop_gt;
        vector<float> percent(HaplotypePath::dbg_mapstop.size());
        for (size_t i = 0; i < HaplotypePath::dbg_mapstop.size(); i++)
            percent[i] = (HaplotypePath::dbg_mapstop[i]/(float)total)*100;

        float percent_acc = 0.0;
        cout << "PBWT mappings stopped after (index : count):" << endl;
        for (size_t i = 0; i < HaplotypePath::dbg_mapstop.size(); i++) {
            percent_acc += percent[i];
            cout << "  " << i << " : " << HaplotypePath::dbg_mapstop[i] << " (" << percent[i] << "%, acc: " << percent_acc << "%)" << endl;
        }
        cout << "  >=" << HaplotypePath::dbg_mapstop.size() << " : " << HaplotypePath::dbg_mapstop_gt << " (" << (HaplotypePath::dbg_mapstop_gt/(double)total)*100 << ", acc: 100%)" << endl;
        cout << "  Total: " << total << endl;
    }
    {
        size_t total = 0;
        for (size_t c : HaplotypePath::dbg_mapone)
            total += c;
        total += HaplotypePath::dbg_mapone_gt;
        vector<float> percent(HaplotypePath::dbg_mapone.size());
        for (size_t i = 0; i < HaplotypePath::dbg_mapone.size(); i++)
            percent[i] = (HaplotypePath::dbg_mapone[i]/(float)total)*100;

        float percent_acc = 0.0;
        cout << "PBWT mappings became 1 at (index : count):" << endl;
        for (size_t i = 0; i < HaplotypePath::dbg_mapone.size(); i++) {
            percent_acc += percent[i];
            cout << "  " << i << " : " << HaplotypePath::dbg_mapone[i] << " (" << percent[i] << "%, acc: " << percent_acc << "%)" << endl;
        }
        cout << "  >=" << HaplotypePath::dbg_mapone.size() << " : " << HaplotypePath::dbg_mapone_gt << " (" << (HaplotypePath::dbg_mapone_gt/(double)total)*100 << ", acc: 100%)" << endl;
        cout << "  Total: " << total << endl;
    }
#endif
}
