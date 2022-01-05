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

#ifdef USE_CUDA_GPU

#include <iomanip>
//// DEBUG
//#include <fstream>
//// __DEBUG

#include "utils.h"
#include "MyMalloc.h"

#include "GPUHandler.h"

using namespace std;
using namespace hybridsys;

GPUHandler::GPUHandler(Hybridsys &hysys_, const VCFData &vcfdata_)
            : hysys(hysys_),
              engines(),
//              bufferFactory(factory_),
              vcfdata(vcfdata_),
              M(vcfdata_.getNSNPs())
{

    int num_gpus = hysys.getGPUs().size();
    for(int i = 0; i < num_gpus; i++)
        engines.emplace_back(hysys.getGPU(i).getIndex());
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
    size_t nt = DEBUG_TARGET_STOP - DEBUG_TARGET_START + 1;

#else
    size_t nt = vcfdata_.getNTarget();
#endif
    totalTargets = nt;
    targetsRemainingPreProcess = nt;
    targetsRemainingProcessCref = nt;
//    targetsRemainingPcR = nt;
//    targetsRemainingPcPBWT = nt;

}

// return the block size
size_t GPUHandler::initialize(size_t K_, size_t Nrefhapscorr_, size_t NrefhapsAD_) {
    K = K_;
    Kwords = roundToMultiple(K_, UNITWORDS * sizeof(data_type) * 8) / (8 * sizeof(data_type));
    Nrefhapscorr = Nrefhapscorr_;
    NrefhapsAD = NrefhapsAD_;

    for (auto &engine : engines)
        engine.initialize(K, Kwords, M, Nrefhapscorr, NrefhapsAD, vcfdata.getReferenceT_const()[0].getCapacity(), vcfdata.getReferenceT_const()[0].getData());

    gpuBlockSize = engines[0].getBlockSize();
    totalBlocks = divideRounded(totalTargets, gpuBlockSize);

    return gpuBlockSize;
}


void GPUHandler::preProcessGPUOnly(
    tbb::concurrent_bounded_queue<PBWTPhaser::targetQueue_type> &inqueue,
    tbb::concurrent_bounded_queue<PBWTPhaser::cPBWTQueue_type> &outqueue,
    atomic_bool& termination_request,
    int threadIndex
    ) {

    ThreadUtils::setThreadName("eglimp prep GPU");
    cout << "Launched preprocessing on GPU " << threadIndex << "." << endl;

    size_t blocksRem_local = totalBlocks; // not for multiple provide threads! (use global atomic blocksRemLocal then!)
    while(!termination_request && blocksRem_local > 0) {

        PBWTPhaser::targetQueue_type block;
        try {
//            // DEBUG
//            cout << "Prep GPU " << threadIndex << ": Waiting... rem: " << blocksRem_local << endl;

            inqueue.pop(block);
            blocksRem_local--;

//            // DEBUG
//            cout << "Prep GPU " << threadIndex << ": Got #" << (totalBlocks-blocksRem_local) << "/" << totalBlocks << "." << endl;

        } catch (tbb::user_abort &e) {
            if (blocksRem_local > 0)
                cout << "Prep GPU " << threadIndex << ": User abort." << endl;
            break;
        }

//        // DEBUG
//        cout << "Phasing preliminaries." << endl;

        vector<const uint64_t*> tgtData(block->size());
        vector<const uint64_t*> splitData(block->size());
        vector<const uint64_t*> besthapData(block->size());
        vector<uint32_t> crefTOffsets(block->size()+1); // initialized with zeros
        vector<uint32_t*> cPBWTData(block->size());
//        vector<vector<BooleanVector>*> crefTptrs(block->size());
        vector<vector<int>*> count0s(block->size());

#pragma omp parallel for
        for (size_t t = 0; t < block->size(); t++) {
            block->at(t)->prepPhase();
            tgtData[t] = block->at(t)->getTargetGTs().getData().data();
            splitData[t] = block->at(t)->getSplits().getData();
            besthapData[t] = block->at(t)->getBestHapsFlags().getData();
            size_t ncref = 2*(block->at(t)->getNSplits())+1;
            crefTOffsets[t+1] = ncref; // due to parallelization we don't store the accumulation yet, the first element will be left as zero (as the first offset is always zero)
            cPBWTData[t] = (uint32_t*) MyMalloc::malloc(8*ncref*Kwords*sizeof(uint32_t), string("cPBWTData_t")+to_string(t)+"_br"+to_string(blocksRem_local));
//            crefTptrs[t] = new vector<BooleanVector>(ncref, BooleanVector());
//            auto curr_dataT = crefTData[t];
////            cout << "t: " << t << " malloc: " << hex << (size_t)crefTData[t] << dec << endl;
//            for (auto &crefT : *(crefTptrs[t])) {
//                crefT.setData(curr_dataT, Kwords*sizeof(uint64_t), K); // size is set to K already, contents will be filled by the GPU kernel
//                curr_dataT += Kwords;
//            }
            count0s[t] = new vector<int>(ncref);
        }

        // accumulate offsets (GPU needs this for the destination data offset)
        uint32_t curroff = 0;
        for (auto &off : crefTOffsets) {
            curroff += off;
            off = curroff;
        }

//        // DEBUG
//        cout << "Running GPU kernel." << endl

        engines[threadIndex].runCrefPBWTKernel(tgtData, splitData, besthapData, crefTOffsets, cPBWTData, count0s);

//        // DEBUG dump condensed reference
//        for (size_t t = 0; t < block->size(); t++) {
//            ofstream ofs(string("crefgpu")+to_string(block->at(t)->getIdx()));
//            ofs << "K: " << K << " Kwords: " << Kwords << endl;
//            ofs << "Mspl: " << block->at(t)->getNSplits() << " Mpbwt: " << crefTptrs[t]->size() << endl;
//            data_type* dptr = crefTptrs[t]->at(0).getData();
//            for (size_t m = 0; m < crefTptrs[t]->size(); m++) {
//                ofs << m << ":";
//                if (m%2)
//                    ofs << " tgtsite: " << block->at(t)->getSplitSites()[m/2];
//                else
//                    ofs << " inc span: " << (((m ==crefTptrs[t]->size()-1) ? M : block->at(t)->getSplitSites()[m/2]) - (m ? (block->at(t)->getSplitSites()[m/2-1]+1) : 0));
//                ofs << hex;
//                for (size_t k = 0; k < Kwords; k++) {
//                    if (k%8 == 0)
//                        ofs << endl;
//                    ofs << " " << setw(16) << setfill('0') << *dptr++;
//                }
//                ofs << dec << endl;
//            }
//            ofs.close();
//        }
//        // __DEBUG

//        // DEBUG dump PBWT
//        for (size_t t = 0; t < block->size(); t++) {
//            ofstream ofs(string("pbwtgpu")+to_string(block->at(t)->getIdx()));
//            ofs << "K: " << K << " Kpbwt: " << 4*Kwords << " Kpad: " << (4*Kwords - 2*divideRounded((int)K,32)) << endl;
//            size_t Mpbwt = block->at(t)->getNSplits() * 2 + 1;
//            ofs << "Mspl: " << block->at(t)->getNSplits() << " Mpbwt: " << Mpbwt << endl;
//            uint32_t* dptr = cPBWTData[t];
//            for (int fwdbck = 0; fwdbck < 2; fwdbck++) {
//                if (fwdbck == 0)
//                    ofs << "fwd:" << endl;
//                else
//                    ofs << "bck:" << endl;
//                for (size_t midx = 0; midx < Mpbwt; midx++) {
//                    size_t m = fwdbck ? Mpbwt-midx-1 : midx;
//                    ofs << m << ":";
////                    if (m%2)
////                        ofs << " tgtsite: " << block->at(t)->getSplitSites()[m/2];
////                    else
////                        ofs << " inc span: " << (((m == Mpbwt-1) ? M : block->at(t)->getSplitSites()[m/2]) - (m ? (block->at(t)->getSplitSites()[m/2-1]+1) : 0));
//                    int count0 = (*(count0s[t]))[m];
//                    int count0x = dptr[2*divideRounded((int)K,32)-1];
//                    ofs << " count0: " << count0 << (count0 == count0x ? " ok" : " XXX");
//                    for (size_t k = 0; k < 2*Kwords; k++) {
//                        if (k%8 == 0)
//                            ofs << endl;
//                        ofs << " " << hex << setw(8) << setfill('0') << *dptr++; // data
//                        ofs << " " << dec << setw(5) << *dptr++; // offset
//                    }
//                    ofs << endl;
//                }
//            }
//            ofs.close();
//        }
//        // __DEBUG

        // insert all targets together with their condensed references into the output queue
        for (size_t t = 0; t < block->size(); t++) {
            PBWTPhaser::cPBWTQueue_type tgt;
            tgt.target = block->at(t);
            tgt.cpbwt = cPBWTData[t];
            tgt.count0 = count0s[t];
            outqueue.push(move(tgt));
        }

    } // end while()

    engines[threadIndex].free();

    if (blocksRem_local > 0)
        cout << "Prep GPU terminated." << endl;
    else
        cout << "Prep GPU finished." << endl;

}

void GPUHandler::processCondensedRef(
    tbb::concurrent_bounded_queue<PBWTPhaser::cPBWTQueue_type> &inqueue,
    tbb::concurrent_bounded_queue<PBWTPhaser::confidence_type> &outqueue,
    atomic_bool& termination_request,
    int threadIndex
    ) {

    ThreadUtils::setThreadName("eglimp cref GPU");
    cout << "Launched processing crefs GPU " << threadIndex << "." << endl;

    while (!termination_request && targetsRemainingProcessCref > 0) {
        PBWTPhaser::cPBWTQueue_type target;
        try {
//            // DEBUG
//            cout << "cref GPU " << threadIndex << ": Waiting... rem: " << targetsRemainingProcessCref << endl;

            inqueue.pop(target);

//            // DEBUG
//            cout << "cref GPU " << threadIndex << ": Got #" << (totalTargets-targetsRemainingProcessCref+1) << "/" << totalTargets << "." << endl;
        } catch (tbb::user_abort &e) {
            if (targetsRemainingProcessCref > 0)
                cout << "cref GPU " << threadIndex << ": User abort." << endl;
            break;
        }

        Target* tgt = target.target;
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT
        cout << "Phasing target " << tgt->getIdx() << "/" << totalTargets << endl;
#endif
        PBWT pbwt(target.cpbwt, target.count0, K, 4*Kwords, 2*target.target->getNSplits()+1, target.target->getIdx());
        tgt->phaseExtPBWT(pbwt);
        PBWTPhaser::confidence_type tconf;
        tconf.id = tgt->getIdx();
        tconf.totalconf = tgt->getTotalConfidence();
        tconf.ncalls = tgt->getNCallSites();
        outqueue.push(tconf);

//        // cleanup memory used for reference and incon data for this target
//        // the memory is shared and starts with the first incon, so freeing the data pointer for the first incon is enough
//        free((*target.cpbwt)[0]);
        // NOTE: target.cpbwt and target.count0 are cleared by the destructor of PBWT

        delete target.target;
//        delete target.cpbwt;
//        delete target.count0;
        targetsRemainingProcessCref--;
    }

    if (targetsRemainingProcessCref > 0)
        cout << "Processing crefs from GPU terminated." << endl;
    else
        cout << "Processing crefs from GPU finished." << endl;

}

#endif // USE_CUDA_GPU
