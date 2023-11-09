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

#ifdef USE_AD_FPGA

#include <array>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <utility>
#include <iostream>
#include <iomanip>
#include <exception>

#include "hybridsys/Hybridsys.h"
#include "hybridsys/ThreadUtils.h"

#include "Datatypes.h"
#include "FPGAConfigurationEagleImp.h"
#include "utils.h"
#include "PBWTPhaser.h"
#include "GPUHandler.h"
#include "StatusFile.h"
#include "MyMalloc.h"
#include "Stopwatch.h"

#include "FPGAHandler.h"

using namespace std;
using namespace hybridsys;

FPGAHandler::FPGAHandler(
        Hybridsys &hysys_,
        chrono::milliseconds timeout_,
        BufferFactory<FPGABuffer> &fpgafactory_wr,
        BufferFactory<FPGABuffer> &fpgafactory_rd,
        BufferFactory<CUDABuffer> &gpufactory,
        const VCFData &vcfdata_,
        size_t maxrefincsites_,
        unsigned numthreads_,
        bool fpgaOnly_,
        bool debug_)
    : hysys(hysys_), timeout(timeout_),
      bufferFactoryFPGA_write(fpgafactory_wr),
      bufferFactoryFPGA_read(fpgafactory_rd),
      bufferFactoryGPU(gpufactory),
      vcfdata(vcfdata_), maxpbwtsites(maxrefincsites_),
      numthreads(numthreads_),
      fpgaOnly(fpgaOnly_), debug(debug_)
{
    hysys.getFPGA(0).createThreadHandle();
    fpgaconf.parse(hysys.getFPGA(0).getConfigurationData());
    hysys.getFPGA(0).destroyThreadHandle();

    blocksize = fpgaconf.getNumPipelines();

    // initialize DMA channels
    for(FPGA &f : hysys.getFPGAs()) {
        f.declareDMAChannel(0, FPGA::DMAChannelType::AxiStreamMaster);
        if (debug)
            cout << "FPGA: " << f.getIndex() << " Channel 0: Master" << endl;
        f.declareDMAChannel(1, FPGA::DMAChannelType::AxiStreamSlave);
        if (debug)
            cout << "FPGA: " << f.getIndex() << " Channel 1: Slave" << endl;

        // we take the lock later
//        Stopwatch swlockfpga("Wait lock (FPGA)");
//        cout << "Waiting for lock on FPGA #" << f.getIndex() << " (sn: " << f.getSerialNumber() << ")... " << flush;
//        f.lock();
//        cout << "got it." << endl;
//        swlockfpga.stop();
    }
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
    size_t nt = DEBUG_TARGET_STOP - DEBUG_TARGET_START + 1;

#else
    size_t nt = vcfdata.getNTarget();
#endif
    totalTargets = nt;
    targetsRemaining_in = nt;
    targetsRemaining_out = nt;
    totalBlocks = nt / blocksize + (nt % blocksize ? 1 : 0);
    fpga_blocks_sent = 0;
    last_block_flag = false;

    // Check if buffers will be large enough

    // buffer must be at least large enough to provide the constants
    constrequired_bufsize = num_constants * sizeof(uint32_t) + fpgaconf.getNumPipelines() * 256/8;
    if(bufferFactoryFPGA_write.getBufferSize() < constrequired_bufsize) {
        StatusFile::addError("Buffer is not large enough to hold the FPGA initialization constants (need at least " + std::to_string(constrequired_bufsize) + " bytes)");
        exit(EXIT_FAILURE);
    }

    // buffer must also be large enough to keep at least one target
    tgtsize = vcfdata.getTargets()[0].getData().size() * sizeof(GenotypeVector::data_type);
    if (bufferFactoryFPGA_write.getBufferSize() < tgtsize) {
        StatusFile::addError("Buffer is not large enough to carry at least one target genotype. Need at least " + std::to_string(tgtsize) + " bytes.");
        exit(EXIT_FAILURE);
    }

//    dbgbufs[0] = (char*) malloc(bufferFactory.getBufferSize());
//    dbgbufs[1] = (char*) malloc(bufferFactory.getBufferSize());
//    dbgbufs[2] = (char*) malloc(bufferFactory.getBufferSize());
}

void FPGAHandler::initIteration(size_t K_, size_t iter_) {
    K = K_;
    iter = iter_;
    targetsRemaining_in = totalTargets;
    targetsRemaining_out = totalTargets;

    Stopwatch swlockfpga("Wait lock (FPGA)");
    for(FPGA &f : hysys.getFPGAs()) { // this should be only one, otherwise this is a perfect deadlock situation...
        cout << "Waiting for lock on FPGA #" << f.getIndex() << " (sn: " << f.getSerialNumber() << ")... " << flush;
        f.lock();
        cout << "got it." << endl;
    }
    swlockfpga.stop();
}

void FPGAHandler::prepareReferenceBuffers() {
    // prepare buffers with ref data
    const size_t singlesitesize = roundToMultiple(iter > 1 ? vcfdata.getNReferenceHapsMax() : vcfdata.getNReferenceHaps(), (size_t)512) / 8;
    const size_t fullrefsize   = vcfdata.getNSNPs() * singlesitesize;
    const size_t bufsize = bufferFactoryFPGA_write.getBufferSize();
    const size_t numbufs = divideRounded(fullrefsize, bufsize);

    if (debug)
        cout << "singlesitesize: " << singlesitesize << " fullrefsize: " << fullrefsize << " bufsize: " << bufsize << " numbufs: " << numbufs << endl;

    refbuffers.clear();
    refbuffers.reserve(2*numbufs); // forward and reverse run
    for (int reverse = 0; reverse < 2; reverse++) {
        for (size_t bidx = 0; bidx < numbufs; bidx++) {
            refbuffers.emplace_back(bufsize);
            refbuffers.back().setContentLength(bufsize);
        }
        if (fullrefsize % bufsize) // last buffer is not full
            refbuffers.back().setContentLength(fullrefsize % bufsize); // correct content length
    }

    size_t currbuf = 0;
    char *curr_bufptr = refbuffers[currbuf].getData();
    size_t buf_remsize = refbuffers[currbuf].getContentLength();

    for (int reverse = 0; reverse < 2; reverse++) { // fill forward and then reverse buffers
        for (size_t m = 0; m < vcfdata.getNSNPs(); m++) {
            size_t copysize = min(buf_remsize, singlesitesize);
            const char* srcptr = reinterpret_cast<const char*>(vcfdata.getReferenceT_const()[reverse ? vcfdata.getNSNPs()-m-1 : m].getData());

            if (copysize > 0) { // memcpy documentation says, src and dest have to be valid even when count is 0
                // add reference to buffer
                memcpy(curr_bufptr, srcptr, copysize);
                curr_bufptr += copysize;
                buf_remsize -= copysize;
            }

            if (copysize != singlesitesize) { // partially copied reference only
                // buffer is full now -> take next buffer
                currbuf++;
                // reset buffer pointer
                curr_bufptr = refbuffers[currbuf].getData();
                buf_remsize = refbuffers[currbuf].getContentLength();
                // have not incremented src ptr yet
                srcptr += copysize;
                // copy rest of reference
                memcpy(curr_bufptr, srcptr, singlesitesize-copysize);
                curr_bufptr += singlesitesize-copysize;
                buf_remsize -= singlesitesize-copysize;
            }
        }
    }
}


void FPGAHandler::initTargets(const vector<Target*> &targets, bool lastblock, int threadIndex) {

    // 1. initialize global constants

    // prepare transmission buffer
    if (debug)
        cout << "Provide FPGA " << threadIndex << ": Waiting for buffer." << endl;
    auto buffer = bufferFactoryFPGA_write.get();
    if (debug)
        cout << "Provide FPGA " << threadIndex << ": Got buffer: " << hex << setw(16) << setfill('0') << (size_t)(buffer->getData()) << dec << endl;

    // check if current target block fits into FPGA memory
    size_t accpbwtsites = 0;
    for (const auto &tgt : targets)
        accpbwtsites += 2*tgt->getNSplits()+1;
    // TODO this is not good behaviour! Better reduce the number of targets in a block in advance to prevent this error from happening!
    if (accpbwtsites > maxpbwtsites) {
        StatusFile::addError("FPGA memory is too low to hold " + to_string((accpbwtsites-targets.size())/2) + " split sites (max: " + to_string((maxpbwtsites-targets.size())/2) + ")");
        exit(EXIT_FAILURE);
    }

    uint32_t *constants = reinterpret_cast<uint32_t*>(buffer->getData());

    // sync word
    constants[0] = 0xDEADBEEFU;
    constants[1] = 0xDEADBEEFU;
    constants[2] = 0xDEADBEEFU;
    constants[3] = 0xDEADBEEFU;
    constants[4] = 0xDEADBEEFU;
    constants[5] = 0xDEADBEEFU;
    constants[6] = 0xDEADBEEFU;
    constants[7] = 0xDEADBEEFU;

    // general target information:
    // # words per target -1 (last target word idx for each target) with split site information (1bit) (target genotype information is 2x this size, calculated on FPGA)
    constants[8] = divideRounded(vcfdata.getNSNPs(), (size_t)512)-1;
    // K PBWT unit words
    constants[9] = K;
    constants[10] = divideRounded(K, (size_t)512); // one unit contains 16x32 bit data

    // general reference information:
    // # total reference words -1 (last ref word idx)
    constants[11] = divideRounded(iter > 1 ? vcfdata.getNReferenceHapsMax() : vcfdata.getNReferenceHaps(), (size_t) 512) * vcfdata.getNSNPs() -1;
    // # last haplotype word per site (512bit words per site -1)
    constants[12] = divideRounded(iter > 1 ? vcfdata.getNReferenceHapsMax() : vcfdata.getNReferenceHaps(), (size_t) 512) -1;

    // general information:
    // # sites (M) -1 (last site)
    constants[13] = vcfdata.getNSNPs()-1;
    // buffer size
    constants[14] = bufferFactoryFPGA_read.getBufferSize() / (256/8); // buffer size of the result buffer in 256bit PCIe words
    // mark the last charge
    constants[15] = lastblock ? 0x80000000U : 0;

    // now the required counts for each target
    uint32_t *tgtcounts = &constants[num_constants];

    for (size_t t = 0; t < targets.size(); t++) {
        // 8 uint32_t's per bus initialization word (256 bit):
        // 1st: no. of split sites
        // 2nd: no. of ref/incon segment words (RAM words required for one condensed reference)
        // others: unused (5 to 7 are not even provided to the initialization pipeline)
        tgtcounts[8*t]   = targets[t]->getNSplits();
        tgtcounts[8*t+1] = (2*targets[t]->getNSplits()+1)*divideRounded(K, (size_t)512) * 2; // required RAM words for one PBWT (2*#splits+1 refincon sites) each K/512 bit PBWT RAM words -> times 2 for fwd/bck
        // not necessarily required:
        for (unsigned x = 2; x < 8; x++)
            tgtcounts[8*t+x] = 0; // unused
    }
    // dummy targets (if number of targets is less than the blocksize)
    for (unsigned t = targets.size(); t < blocksize; t++) { // 4 uint32_t's per bus initialization word
        tgtcounts[8*t]   = 0; // no split sites for dummy targets
        tgtcounts[8*t+1] = divideRounded(K, (size_t)512) * 2; // but each dummy target produces a PBWT with one site with K/512 bit PBWT RAM words that have to be recognized -> times 2 for fwd/bck
        // not necessarily required:
        for (unsigned x = 2; x < 8; x++)
            tgtcounts[8*t+x] = 0; // unused
    }

    if (debug) {
        cout << "Provide FPGA " << threadIndex << ": FPGA initialization of targets with constants: " << constrequired_bufsize << " bytes." << endl;
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
        cout << hex;
        uint32_t* bufptr = constants;
        for (size_t x = 0; x < constrequired_bufsize/4; x++) {
            cout << " " << setw(8) << setfill('0') << *bufptr++ << " ";
            if (x % 8 == 7)
                cout << endl;
        }
        cout << dec << endl;
#endif
    }

    // NOTE: constants will be transmitted together with first target charge
    char *curr_bufptr = buffer->getData() + constrequired_bufsize;
    size_t buf_remsize = bufferFactoryFPGA_write.getBufferSize() - constrequired_bufsize;

    { // target data
        if (debug) {
            cout << "Provide FPGA " << threadIndex << ": Target data: single: " << tgtsize << " total: " << (blocksize*tgtsize) << endl;
        }
        for (size_t tidx = 0; tidx < blocksize; tidx++) {
            if (buf_remsize < tgtsize) {
                // not enough space in buffer for next target -> send to FPGA
                if (debug)
                    cout << "Provide FPGA " << threadIndex << ": (Target) data to FPGA: " << bufferFactoryFPGA_write.getBufferSize() - buf_remsize << " bytes." << endl;

                hysys.getFPGA(threadIndex).writeDMA(*buffer, 0, timeout, bufferFactoryFPGA_write.getBufferSize() - buf_remsize);
                // reset buffer pointer
                curr_bufptr = buffer->getData();
                buf_remsize = bufferFactoryFPGA_write.getBufferSize();
            }
            if (tidx < targets.size()) {
                // add target to buffer
                size_t curr_target = targets[tidx]->getIdx();
                memcpy(curr_bufptr, reinterpret_cast<const char*>(vcfdata.getTargets()[curr_target].getData().data()), tgtsize);
            } // else leave buffer content with garbage for non-existing targets
            curr_bufptr += tgtsize;
            buf_remsize -= tgtsize;
        }
    }
    // buffer still contains data to be sent!

    { // split data
        size_t splitsize = targets[0]->getSplits().getDataSize();
        if (debug) {
            cout << "Provide FPGA " << threadIndex << ": Split data: single: " << splitsize << " total: " << (blocksize*splitsize) << endl;
        }
        for (size_t tidx = 0; tidx < blocksize; tidx++) {
            if (buf_remsize < splitsize) {
                // not enough space in buffer for next target split sites -> send to FPGA
                // DEBUG
                if (debug)
                    cout << "Provide FPGA " << threadIndex << ": (Split) data to FPGA: " << bufferFactoryFPGA_write.getBufferSize() - buf_remsize << " bytes." << endl;

                hysys.getFPGA(threadIndex).writeDMA(*buffer, 0, timeout, bufferFactoryFPGA_write.getBufferSize() - buf_remsize);
                // reset buffer pointer
                curr_bufptr = buffer->getData();
                buf_remsize = bufferFactoryFPGA_write.getBufferSize();
            }
            if (tidx < targets.size()) {
                // add target split sites to buffer
                memcpy(curr_bufptr, reinterpret_cast<const char*>(targets[tidx]->getSplits().getData()), splitsize);
            } else {
                // fill buffer with zeros for non-existing targets
                memset(curr_bufptr, 0, splitsize);
            }
            curr_bufptr += splitsize;
            buf_remsize -= splitsize;
        }
    }
    // buffer still contains data to be sent!

    { // best haps data
        size_t bhsize = targets[0]->getBestHapsFlags().getDataSize();
        if (debug) {
            cout << "Provide FPGA " << threadIndex << ": BH data: single: " << bhsize << " total: " << (blocksize*bhsize) << endl;
        }
        for (size_t tidx = 0; tidx < blocksize; tidx++) {
            if (buf_remsize < bhsize) {
                // not enough space in buffer for next target split sites -> send to FPGA
                // DEBUG
                if (debug)
                    cout << "Provide FPGA " << threadIndex << ": (Best haps) data to FPGA: " << bufferFactoryFPGA_write.getBufferSize() - buf_remsize << " bytes." << endl;

                hysys.getFPGA(threadIndex).writeDMA(*buffer, 0, timeout, bufferFactoryFPGA_write.getBufferSize() - buf_remsize);
                // reset buffer pointer
                curr_bufptr = buffer->getData();
                buf_remsize = bufferFactoryFPGA_write.getBufferSize();
            }
            if (tidx < targets.size()) {
                // add target best haps to buffer
                memcpy(curr_bufptr, reinterpret_cast<const char*>(targets[tidx]->getBestHapsFlags().getData()), bhsize);
                curr_bufptr += bhsize;
            } else {
                // fill buffer with K ones and rest zeros for non-existing targets
                size_t onebytes = K/8;
                if (onebytes) { // bytes completely filled with ones
                    memset(curr_bufptr, 0xffu, onebytes);
                    curr_bufptr += onebytes;
                }
                if (K%8) { // the remaining ones
                    memset(curr_bufptr, (1u<<(K%8))-1u, 1);
                    curr_bufptr++;
                    onebytes++;
                }
                if (bhsize-onebytes) { // remaining zeros
                    memset(curr_bufptr, 0x00u, bhsize-onebytes);
                    curr_bufptr += bhsize-onebytes;
                }
            }
            buf_remsize -= bhsize;
        }
    }
    // buffer still contains data to be sent!

    // send remaining data in buffer
    // DEBUG
    if (debug)
        cout << "Provide FPGA " << threadIndex << ": Remaining const/tgt/split data to FPGA: " << bufferFactoryFPGA_write.getBufferSize() - buf_remsize << " bytes." << endl;

//    if (bufferFactoryFPGA.getBufferSize() - buf_remsize > 0) { // DEBUG: not necessary without debugging, buffer will always contain data then
        FPGA::Result result = hysys.getFPGA(threadIndex).writeDMA(*buffer, 0, timeout, bufferFactoryFPGA_write.getBufferSize() - buf_remsize);
        if (result == FPGA::Result::Timeout) {
            cout << endl;
            StatusFile::addError("Init FPGA timeout.");
            exit(EXIT_FAILURE);
        }
        else if (result == FPGA::Result::Cancelled) {
            cout << endl;
            StatusFile::addError("Init FPGA cancelled.");
            exit(EXIT_FAILURE);
        }
//    }

    // initialization finished, now the FPGA awaits the references to be streamed

}

void FPGAHandler::processBlock(const vector<Target*> &target_block, bool lastblock, int threadIndex) {
    // initialize targets on FPGA
    initTargets(target_block, lastblock, threadIndex);

    // stream prepared references to the FPGA
    if (debug)
        cout << "Provide FPGA " << threadIndex << ": Reference data to FPGA..." << endl;
    for (const auto &rb : refbuffers) {
        FPGA::Result result;
        result = hysys.getFPGA(threadIndex).writeDMA(rb, 0, timeout, rb.getContentLength());
        if (result == FPGA::Result::Timeout) {
            cout << endl;
            StatusFile::addError("Provide FPGA timeout.");
            exit(EXIT_FAILURE);
        } else if (result == FPGA::Result::Cancelled) {
            cout << endl;
            StatusFile::addError("Provide FPGA cancelled.");
            exit(EXIT_FAILURE);
        }
    }
}

// thread function that prepares the targets to be processed on the FPGA
void FPGAHandler::preparePhasing(
        tbb::concurrent_bounded_queue<Target*> &inqueue,
        tbb::concurrent_bounded_queue<Target*> &outqueue,
        atomic_bool& termination_request,
        int threadIndex
        ) {

    ThreadUtils::setThreadName("eimp prep FPGA");
    if (debug)
        cout << "Launched PrepPhasing FPGA " << threadIndex << "." << endl;

    while (!termination_request && targetsRemaining_in > 0) {

        try {
            Target* tgtptr;
            inqueue.pop(tgtptr);
            targetsRemaining_in--;

            if (debug)
                cout << "PrepPhasing " << threadIndex << ": Preparing target #" << (totalTargets-targetsRemaining_in) << endl; // probable side effect with counter (if another thread has decremented it meanwhile), but debug output anyway

            // preparation + push
            tgtptr->prepPhase();
            outqueue.push(tgtptr);

        } catch (tbb::user_abort &e) {
            if (debug && targetsRemaining_in > 0)
                cerr << "PrepPhasing " << threadIndex << ": User abort." << endl;
            break;
        }

    }

}

// thread function that initializes the targets and streams the references
void FPGAHandler::provideReferences(
        tbb::concurrent_bounded_queue<Target*> &inqueue,
        tbb::concurrent_bounded_queue<PBWTPhaser::targetQueue_type> &outqueue,
        atomic_bool& termination_request,
        int threadIndex
        ) {

    ThreadUtils::setThreadName("eimp prov FPGA");
    if (debug)
        cout << "Launched Provide FPGA " << threadIndex << "." << endl;

    size_t blocksRem_local = totalBlocks; // not for multiple provide threads! (use global atomic blocksRemLocal then!)
    size_t tgtsRem_local = totalTargets;
    size_t curr_filled = 0;

    vector<Target*> *target_block = new vector<Target*>;
    while (tgtsRem_local > 0) { // we need to handle a termination request inside the loop and break

        Target* tgtptr;

        try {
            if (!termination_request)
                inqueue.pop(tgtptr);
            else
                throw tbb::user_abort();

        } catch (tbb::user_abort &e) {
            if (debug && tgtsRem_local > 0)
                cerr << "Provide FPGA " << threadIndex << ": User abort." << endl;
            break; // leave the loop with a termination request
        }

        // need to push the current data to the queue,
        // in the case of a termination request, we know at least one datum is in the block!
        if (curr_filled == blocksize) {
            blocksRem_local--;
            if (debug)
                cout << "Provide FPGA " << threadIndex << ": Got block #" << (totalBlocks-blocksRem_local) << "/" << totalBlocks << " (tgtsrem: " << tgtsRem_local << ")." << endl;

            processBlock(*target_block, false, threadIndex); // lastblock=false: this can never be the last block as there is already one target that needs to be pushed into a new block at this point

            fpga_blocks_sent++; // set this before actually sending the block to prevent interprocess side effects
            outqueue.push(PBWTPhaser::targetQueue_type(target_block)); // make a shared pointer from target_block and add to queue

            if (debug) {
                cout << "Provide FPGA " << threadIndex << ": Pushed #" << (totalBlocks - blocksRem_local) << "/" << totalBlocks <<  " (size: " << target_block->size() << ")." << endl;
            }

            target_block = new vector<Target*>; // assign new block
            curr_filled = 0;
        }

        target_block->push_back(tgtptr);
        curr_filled++;
        tgtsRem_local--;

    }

    // send the last block to the FPGA
    // (even with a termination request, we know that there's at least one target in the block)
    blocksRem_local--;
    if (debug)
        cout << "Provide FPGA " << threadIndex << ": Got last block #" << (totalBlocks-blocksRem_local) << "/" << totalBlocks << " (tgtsrem: " << tgtsRem_local << ")." << endl;

    processBlock(*target_block, true, threadIndex);

    fpga_blocks_sent++;
    last_block_flag = true;
    outqueue.push(PBWTPhaser::targetQueue_type(target_block)); // make a shared pointer from target_block and add to queue

    if (debug) {
        cout << "Provide FPGA " << threadIndex << ": Pushed last block #" << (totalBlocks - blocksRem_local) << "/" << totalBlocks <<  " (size: " << target_block->size() << ")." << endl;
    }

    if (debug) {
        if (tgtsRem_local > 0)
            cout << "Provide FPGA " << threadIndex << " terminated." << endl;
        else
            cout << "Provide FPGA " << threadIndex << " finished." << endl;
    }

}


//void FPGAHandler::readFPGA(
//             atomic_bool& termination_request,
//             int threadIndex
//             ) {
//
//
//    FPGA &fpga = hysys.getFPGA(threadIndex);
//
//    ThreadUtils::setThreadName("eimp read FPGA");
//    if (debug)
//        cout << "Launched Read FPGA " << threadIndex << "." << endl;
//
//    // for debug
//    const size_t buffer512words = (bufferFactoryFPGA_read.getBufferSize()*8)/512; // number of 512bit words in buffer
//    size_t buffersread = 0;
//
//    while(!termination_request) {
//
//        if (debug)
//            cout << "Read FPGA " << threadIndex << ": Waiting for buffer..." << endl;
//
//        shared_ptr<FPGABuffer> b = bufferFactoryFPGA_read.get();
//
//        if (debug) {
//            cout << "Read FPGA " << threadIndex << ": Reading from FPGA...(buffer " << buffersread << ") " << endl;
//            buffersread++;
//        }
//
//        // read the buffer from FPGA
//        FPGA::Result result = fpga.readDMA(*b, 1, timeout); // don't let this run into a timeout! data loss may occur!
//
//        switch(result) {
//        case FPGA::Result::Success:
//            {
//                if (debug)
//                    cout << "Read FPGA " << threadIndex << ": Received " << buffer512words << " in buffer." << endl;
//                fpgaReadQueue.push(b);
//            }
//            break;
//        case FPGA::Result::Timeout:
//            cout << endl;
//            StatusFile::addError("Read FPGA timeout.");
//            exit(EXIT_FAILURE);
//            break;
//        case FPGA::Result::Cancelled: // normal way to stop this thread TODO how to make this nicer?!
//            if (debug)
//                cout << "Read FPGA cancelled." << endl;
//            break;
//        }
//
//    }
//
//    if (debug)
//        cout << "Finished Read FPGA " << threadIndex << "." << endl;
//}


void FPGAHandler::preProcessFPGAOnly(tbb::concurrent_bounded_queue<PBWTPhaser::targetQueue_type> &inqueue,
             tbb::concurrent_bounded_queue<PBWTPhaser::cPBWTQueue_type> &outqueue,
             atomic_bool& termination_request,
             int threadIndex
             ) {

//    // DEBUG
//    size_t dbgcopycount = 0;
//    size_t dbgnocopycount = 0;

    FPGA &fpga = hysys.getFPGA(threadIndex);

    ThreadUtils::setThreadName("eimp prep FPGA");
    if (debug)
        cout << "Launched PreProc FPGA " << threadIndex << "." << endl;

//    // spawn separate reader thread
//    atomic_bool reader_term_request;
//    reader_term_request = false;
//    thread readthr(&FPGAHandler::readFPGA, this, std::ref(reader_term_request), threadIndex);

    const size_t site512words = divideRounded(K, (size_t)512); // how many 512 bit words per site?
    size_t blocks_processed = 0; // not for multiple threads! (use global atomic then!)
    const size_t buffer512words = (bufferFactoryFPGA_read.getBufferSize()*8)/512; // number of 512bit words in buffer
    size_t buffer_rem512words = 0; // unprocessed 512bit words left in buffer
    shared_ptr<FPGABuffer> b;
    uint32_t* data = NULL;
    size_t buffersread = 0; // for debug

    while(!last_block_flag || blocks_processed < fpga_blocks_sent) { // this ensures that all data from blocks sent to the FPGA will be fetched (especially in the case of a user termination)

        PBWTPhaser::targetQueue_type block;
        try {
            if (debug)
                cout << "PreProc FPGA " << threadIndex << ": Waiting... rem: " << (totalBlocks - blocks_processed) << endl;

            inqueue.pop(block);
            blocks_processed++;

            if (debug)
                cout << "PreProc FPGA " << threadIndex << ": Got #" << blocks_processed << "/" << totalBlocks << "." << endl;

        } catch (tbb::user_abort &e) { // should not occur since the queue should not be allowed to be aborted!!!
            if (debug && blocks_processed < totalBlocks)
                cerr << "PreProc FPGA " << threadIndex << ": User abort." << endl;
            break;
        }

        // the complete condensed PBWT data for this block
        vector<uint32_t*> allpbwts(block->size(), NULL);
//        // reserve for each target the required amount of space for the PBWT according to the number of split sites
//        for (size_t t = 0; t < block->size(); t++) {
//            size_t capacity = site512words*512/8; //roundToMultiple(sitesize, UNITWORDS * sizeof(BooleanVector::data_type) * 8) / 8;
//            size_t nsplits = (*block)[t]->getNSplits();
//            uint32_t *pbwtdata = (uint32_t*) MyMalloc::malloc((2*nsplits+1) * capacity * 2, string("pbwtdata_t")+to_string(t)+"_b"+to_string(blocks_processed)); // space for the fwd and reverse PBWT for this target
//            allpbwts[t] = pbwtdata; //new PBWTPhaser::refincs_type(2*nsplits+1, BooleanVector(refincdata, capacity));
//        }
        // add dummy targets to total number of words, if necessary (should be the case only for last buffer)
        size_t dummytgt512words = (blocksize - block->size()) * site512words * 2; // dummies also produce a fwd and a bck PBWT

        // to keep track of where we are:
        size_t curr_tgt = 0;  // curr target (relative to block)
        uint32_t *curr_data = allpbwts[curr_tgt]; // start of PBWT data block for this target
        size_t tgt_512words = (2*(*block)[curr_tgt]->getNSplits()+1) * site512words * 2; // including fwd and bck PBWT
        size_t tgt_rem512words = tgt_512words;

        // process all required 512 bit words for this block
#if defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET
        size_t w = 0;
#endif
        while (tgt_rem512words > 0) {

            if (buffer_rem512words == 0) { // load next buffer from queue
                if (debug)
                    cout << "PreProc FPGA " << threadIndex << ": Waiting for buffer..." << endl;

//                fpgaReadQueue.pop(b);

                b = bufferFactoryFPGA_read.get();

                if (debug)
                    cout << "PreProc FPGA " << threadIndex << ": Reading from FPGA...(buffer " << buffersread << ") " << endl;

                // read the buffer from FPGA
                FPGA::Result result = fpga.readDMA(*b, 1, timeout); // don't let this run into a timeout! data loss may occur!

                switch(result) {
                case FPGA::Result::Success:
                    {
                        if (debug) {
                            cout << "PreProc FPGA " << threadIndex << ": Received " << buffer512words << " in buffer." << endl;
                            buffersread++;
                        }
//                        fpgaReadQueue.push(b);
                    }
                    break;
                case FPGA::Result::Timeout:
                    cout << endl;
                    StatusFile::addError("PreProc FPGA read timeout.");
                    exit(EXIT_FAILURE);
                    break;
                case FPGA::Result::Cancelled:
                    cout << endl;
                    StatusFile::addError("PreProc FPGA read cancelled.");
                    exit(EXIT_FAILURE);
                    break;
                }

//                if (debug) {
//                    cout << "PreProc FPGA " << threadIndex << ": Got buffer " << buffersread << endl;
//                    buffersread++;
//                }

                data = reinterpret_cast<uint32_t*>(b->getData());
                buffer_rem512words = buffer512words;
            } // end if buffer_rem512words == 0

            // data copying from buffer:
            // process as much data as possible for this target
            size_t copy512words = min(buffer_rem512words, tgt_rem512words);

            bool tgtcomplete = false;
            if (tgt_512words == tgt_rem512words && copy512words == tgt_rem512words && curr_tgt < block->size()) {
                // complete target is in the buffer -> just set the pointers to the data and skip copying
                tgtcomplete = true;
                allpbwts[curr_tgt] = data;
            }

            tgt_rem512words -= copy512words;
            buffer_rem512words -= copy512words;

#if defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET
            if(debug) {
                cout << "Tgt: " << (curr_tgt < block->size() ? to_string(curr_tgt).c_str() : "pad!")
                        << "\nrem: " << (tgt_rem512words+copy512words)
//                        << "\nrem w. pad: " << (tgt_rem512words+copy512words+tgt_rem512words_padding)
//                        << "\ntot rem: " << (total512words - w)
                        << "\nbuf: " << buffer_rem512words
                        << "\nnow copying: " << copy512words << " (" << copy512words*8*sizeof(uint64_t) << " bytes from "
                        << hex << setw(16) << setfill('0') << (uint64_t) data << " to "
                        << setw(16) << setfill('0') << (uint64_t) curr_data << ")" << dec << endl;
            }
#endif

            // copy
            if (!tgtcomplete && copy512words > 0 && curr_tgt < block->size() && !termination_request) { // only copy really required data
                // check if space was already allocated for the destination.
                // if not, allocate!
                if (!curr_data) {
                    size_t capacity = site512words*512/8; //roundToMultiple(sitesize, UNITWORDS * sizeof(BooleanVector::data_type) * 8) / 8;
                    size_t nsplits = (*block)[curr_tgt]->getNSplits();
                    uint32_t *pbwtdata = (uint32_t*) MyMalloc::malloc((2*nsplits+1) * capacity * 2, string("pbwtdata_t")+to_string(curr_tgt)+"_b"+to_string(blocks_processed)); // space for the fwd and reverse PBWT for this target
                    allpbwts[curr_tgt] = pbwtdata; //new PBWTPhaser::refincs_type(2*nsplits+1, BooleanVector(refincdata, capacity));
                    curr_data = pbwtdata;
                }
                memcpy(curr_data, data, copy512words*512/8);
            }

            // increment counters and data ptr
            data += copy512words*512/32;
            curr_data += copy512words*512/32;

#if defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET
            if (debug)
                w += copy512words;
#endif

            // finished a target (or dummies)
            if (tgt_rem512words == 0) {

                // push current target into outqueue
                if (!termination_request && curr_tgt < block->size()) { // only real targets and we don't need to push the data into outqueue if we received a termination request
                    if (debug)
                        cout << "PreProc FPGA " << threadIndex << ": Pushing tgt #" << curr_tgt+1 << "/" << block->size() << " of block #" << blocks_processed << "/" << totalBlocks <<  "." << endl;
                    PBWTPhaser::cPBWTQueue_type tgt;
                    tgt.target = (*block)[curr_tgt];
                    tgt.cpbwt = allpbwts[curr_tgt];
                    if (tgtcomplete)
                        tgt.parent = b; // prevent the buffer from deletion
//                    // DEBUG
//                    if (tgtcomplete)
//                        dbgnocopycount++;
//                    else
//                        dbgcopycount++;
//                    cerr << "CC:" << dbgcopycount << " NC:" << dbgnocopycount << endl;
//                    // __DEBUG
                    outqueue.push(tgt);
                }

                // prepare next target
                curr_tgt++;
                tgtcomplete = false;
                if (curr_tgt < block->size()) { // only real targets
#if defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET
                    if (debug)
                        cout << " Finished target padding. Next target..." << endl;
#endif
                    curr_data = allpbwts[curr_tgt]; // start of PBWT data block for this target (fwd and bck)
                    tgt_512words = (2*(*block)[curr_tgt]->getNSplits()+1) * site512words * 2;
                    tgt_rem512words = tgt_512words;
                } else if (curr_tgt == block->size()){ // padding targets
                    // this only applies for the very last block in this run if we have padding targets,
                    // but is no problem for usual block ends (dummytgt512words is zero then):
                    // we want all remaining dummy words at once now
                    tgt_512words = 0;
                    tgt_rem512words = dummytgt512words;

#if defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET
                    if (debug) {
                        cout << " Finished target padding. Was last target in block." << endl;
//                            if (curr_tgt < blocksize)
//                                cout << " Last block with padding targets! Padding words: " << tgt_rem512words_padding << endl;
                    }
#endif

                } // tgt_rem512words stays zero and we're finished
#if defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET
                    else if (debug)
                        cout << " Finished target padding. Was padding target of block." << endl;
#endif
            }

        } // end while over current block

        // received and processed all buffers for this block

    } // end while over all blocks

//    // stop FPGA reading process
//    reader_term_request = true;
//    try {
//        fpga.cancel(1);
//        readthr.join();
//    } catch (exception& e) {
//        cerr << "Caught an exception while cancelling FPGA read. Continuing anyway..." << endl;
//    }

    if (debug) {
        if (blocks_processed < totalBlocks)
            cout << "PreProc FPGA " << threadIndex << " terminated." << endl;
        else
            cout << "PreProc FPGA " << threadIndex << " finished." << endl;
    }

    // release FPGA lock for other processes
    fpga.unlock();

}


void FPGAHandler::processPBWT(
    tbb::concurrent_bounded_queue<PBWTPhaser::cPBWTQueue_type> &inqueue,
    tbb::concurrent_bounded_queue<PBWTPhaser::confidence_type> &outqueue,
    atomic_bool& termination_request,
    int threadIndex
    ) {

        ThreadUtils::setThreadName("eimp phas FPGA");
        if (debug)
            cout << "Launched Phase FPGA " << threadIndex << "." << endl;

        while (!termination_request && targetsRemaining_out > 0) {
            PBWTPhaser::cPBWTQueue_type target;
            try {
                // DEBUG
                if (debug)
                    cout << "Phase FPGA " << threadIndex << ": Waiting... rem: " << targetsRemaining_out << endl;

                inqueue.pop(target);
                targetsRemaining_out--;

                // DEBUG
                if (debug)
                    cout << "Phase FPGA " << threadIndex << ": Got #" << (totalTargets-targetsRemaining_out) << "/" << totalTargets << "." << endl; // probable side effect regarding the counter, but it's debug output anyway...
            } catch (tbb::user_abort &e) {
                if (debug && targetsRemaining_out > 0)
                    cerr << "Phase FPGA " << threadIndex << ": User abort." << endl;
                break;
            }

            Target* tgt = target.target;
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
            cout << "Phasing target " << tgt->getIdx() << "/" << totalTargets << endl;
#endif

            // number of 32bit words per site
            size_t Kwords = divideRounded(K,(size_t)512)*512/32; // number of 512 bit words to 32bit words
            // number of PBWT sites
            int M = 2*(tgt->getNSplits())+1;

            PBWT pbwt(target.cpbwt, K, Kwords, M, tgt->getIdx());
            if (target.parent.get() == NULL) // no parent -> my memory -> can free it
                MyMalloc::free(target.cpbwt); // has been copied by now
            //else parent exists -> don't free the memory -> is re-used in FPGA-host communication
            tgt->phaseExtPBWT(pbwt);
            PBWTPhaser::confidence_type tconf;
            tconf.id = tgt->getIdx();
            tconf.totalconf = tgt->getTotalConfidence();
            tconf.ncalls = tgt->getNCallSites();
            outqueue.push(tconf);

            // cleanup this target
            delete target.target;
        }

        if (debug) {
            if (targetsRemaining_out > 0)
                cout << "Phase FPGA " << threadIndex << " terminated." << endl;
            else
                cout << "Phase FPGA " << threadIndex << " finished." << endl;
        }
}

#endif // USE_AD_FPGA
