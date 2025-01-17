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

#include <algorithm>
#include <functional>
#include <chrono>
#include <iostream>
#include <fstream>
#include <thread>
#include <omp.h>
#include <csignal>

#include "version.h"
#include "Args.h"
#include "FPGAConfigurationEagleImp.h"
#include "utils.h"
#include "StatusFile.h"
#include "GPUEngine.h"
#include "hybridsys/Hybridsys.h"
#include "hybridsys/ThreadUtils.h"
#include "VCFData.h"
#include "PBWTPhaser.h"
#include "PBWTImputer.h"
#include "Stopwatch.h"

#include "MyMalloc.h"

using namespace std;
using namespace hybridsys;

StatusFile StatusFile::instance;

MyMalloc MyMalloc::mm;

// this is required for a clean shutdown when using FPGAs
static bool term_request = false;
static bool term_wait = false;

/**
 * @brief Signal handler registered for SIGINT reception
 *
 * This function is registered to be called by the operating system if a SIGINT
 * is received (i.e. Ctrl+C). We need to catch this and do a clean process
 * termination when using FPGAs.
 *
 * @param signal The actual signal that has been received
 */
static void signalHandler(int signal) {
    if(signal == SIGINT || signal == SIGTERM) {
        cerr << "Received interrupt. ";
        if(!term_request && term_wait) {
            term_request = true;
            PBWTPhaser::terminate = true;
            cerr << "Trying to shutdown gracefully, please be patient. Hit Ctrl+C again to terminate immediately." << endl;
        } else {
            cerr << "Terminating..." << endl;
            exit(EXIT_FAILURE);
        }
    }
}


int main(int argc, char *argv[]) {

    ThreadUtils::setThreadName("eagleimp");

    cout << "--------------------------------------------------------" << endl;
    cout << "--- EagleImp: Genotype Phasing + Imputation" << endl;
    cout << "---" << endl;
    cout << "--- " << prog_version << endl;
    cout << "--- (compiled on " <<  prog_timestamp << ")" << endl;
    cout << "---" << endl;
    cout << "--- Copyright (C) 2018-2025 by Lars Wienbrandt," << endl;
    cout << "--- Institute of Clinical Molecular Biology, Kiel University" << endl;
    cout << "--- Distributed under the GNU GPLv3 license." << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << endl;

    if (argc > 1) {
        cout << "Provided arguments:" << endl;
        cout << "  ";
        size_t cpos = 0;
        for (int a = 1; a < argc; a++) { // skip first argument (call to binary)
            if (a>1) {
                if (argv[a][0] == '-') {
                    cout << endl << "  ";
                    cpos = 0;
                } else {
                    for (; cpos < 18; cpos++) {
                        cout << ' ';
                    }
                }
            }
            size_t p = string(argv[a]).find_last_of('/');
            if (p == string::npos) {
                cout << argv[a];
                cpos += string(argv[a]).size();
            } else
                cout << "***" << string(argv[a]).substr(p); // mask all paths, do not need to update cpos here
        }
        cout << endl;
    } // else: no arguments provided. will exit in Args.

    Args args = Args::parseArgs(argc, argv);

    { // at the end of this block all destructors for objects inside main were called. Useful for the (debug) usage of MyMalloc.

    // register Ctrl+C handler for proper cleanup on interruption
    signal(SIGINT, signalHandler);

    {
        ifstream f("/proc/self/status");
        string line;
        int tracerPid = 0;
        static const char token[] = "TracerPid:";

        while(getline(f, line)) {
            size_t pos = line.find(token);
            if(pos != string::npos) {
                tracerPid = atoi(line.data() + pos + sizeof(token));
            }
        }
        if (tracerPid > 0)
        	StatusFile::addWarning("Detected attached debugger.");
    }

    // initialize platform
#ifdef USE_AD_FPGA
    vector<unsigned> fpgas = args.get<vector<unsigned>>("fpgas");
#else
    vector<unsigned> fpgas;
#endif
#ifdef USE_CUDA_GPU
    vector<unsigned> gpus = args.get<vector<unsigned>>("gpus");
#else
    vector<unsigned> gpus;
#endif

    bool usefpga, usegpu;
    vector<FPGAConfigurationEagleImp> fpgaconfigs;
    {
        Hybridsys hysys(fpgas, gpus);
        fpgaconfigs = verifyConfiguration(hysys, usefpga, usegpu, args.count("debug"));
    }

    // lockfile
    string lockfile_plain, lockdir;
    if (args.lockfile.size() > 0) { // only if a CPU lock is desired
        if (*(args.lockfile.rbegin()) == '/') {
            StatusFile::addError("Lockfile must be a file, not a path!");
            exit(EXIT_FAILURE);
        }
        if (args.lockfile.find_first_of('/') == string::npos) { // no slash in lockfile -> file relative to working dir
            lockfile_plain = args.lockfile;
            lockdir = ".";
        } else {
            lockdir = args.lockfile.substr(0, args.lockfile.find_last_of('/'));
            if (args.lockfile[0] != '/') // relative path
                lockdir = "./" + lockdir;
            lockfile_plain = args.lockfile.substr(args.lockfile.find_last_of('/')+1);
        }
    }

    if (args.skipPhasing)
        args.iters = 1;
    else if (args.iters == 0)
        args.iters = 1; // TODO automatically set number of iterations -> currently always one

    if (!args.skipImputation) {
        cout << "\nSelected imputation information per genotype:";
        if (args.writeADosage)
            cout << "\n - allele dosages";
        if (args.writeGDosage)
            cout << "\n - genotype dosage";
        if (args.writeProb)
            cout << "\n - genotype posterior probabilities";
        if (!(args.writeADosage || args.writeGDosage || args.writeProb)) {
            cout << " none." << endl;
            StatusFile::addWarning("No imputation info selected. Your imputation output will not contain any dosage or probability information.");
        } else
            cout << endl;
        if (args.setmaxerrors == 0 && (args.overwriteCalls || args.improveCalls)) {
            StatusFile::addWarning("--overwriteCalls and --improveCalls will result in dosages of 0.0 or 1.0 in all called and phased variants,<br>\n"
                    "  as they are used as anchor points for imputation. (Otherwise, they would be set to the phasing dosages.)");
        }
    } else
        args.outputPhased = true; // always write phased targets to output

    // phasing parameters (hidden options)
    Target::cMmaxSplit = args.cMmaxsplit;
//    if (args.min_ibd_a != 0 && args.min_ibd_a < args.min_ibd_b) {
//        cerr << "ERROR! minibdA must be greater or equal to minibdB." << endl;
//        exit(EXIT_FAILURE);
//    }
//    if (args.min_ibd_a == 0 || args.min_ibd_b == 0) { // disable IBD check
//        Target::MinIBDSplitSitesA = 0ull;
//        Target::MinIBDSplitSitesB = 0ull;
//    } else {
//        Target::MinIBDSplitSitesA = args.min_ibd_a;
//        Target::MinIBDSplitSitesB = args.min_ibd_b;
//    }
    if (args.hsize > 128 || args.hsize < args.delta+1) {
        StatusFile::addError("hsize must be <= 128 and >= delta+1.");
        exit(EXIT_FAILURE);
    }
    if (args.hsize_fast > 128 || args.hsize_fast < args.delta_fast+2) {
        StatusFile::addError("hsize_fast must be <= 128 and >= delta_fast+2.");
        exit(EXIT_FAILURE);
    }
    Target::HistorySize = args.hsize;
    Target::HistorySizeFast = args.hsize_fast;
    Target::BeamWidth = args.beamwidth;
    Target::BeamWidthFast = args.beamwidth_fast;
    Target::Delta = args.delta;
    Target::DeltaFast = args.delta_fast;
    if (args.minfixsize > 1.0 || args.minfixsize < 0.0 || args.maxfixsize > 1.0 || args.maxfixsize < 0.0) {
        StatusFile::addError("minfixsize and maxfixsize need to be in the range [0.0,1.0].");
        exit(EXIT_FAILURE);
    }
    Target::minFixSize = args.minfixsize;
    Target::maxFixSize = args.maxfixsize;
    if (args.minfixthr > 1.0 || args.minfixthr < 0.0 || args.maxfixthr > 1.0 || args.maxfixthr < 0.0) {
        StatusFile::addError("minfixthr and maxfixthr need to be in the range [0.0,1.0].");
        exit(EXIT_FAILURE);
    }
    if (args.minfixthr > args.maxfixthr) {
        StatusFile::addError("minfixthr cannot be larger than maxfixthr.");
        exit(EXIT_FAILURE);
    }
    Target::minFixThresh = args.minfixthr;
    Target::maxFixThresh = args.maxfixthr;

    // prepare data (read metadata for automatic distribution in chunks)
    VCFData vcfdata(args, argc, argv, fpgaconfigs);
    usefpga = vcfdata.getUseFPGA(); // VCFData may disable FPGA usage due to FPGA restrictions to the dataset

    uint64_t maxpbwtsites = ~0ull;
    if (usefpga) {
        uint64_t avmem = fpgaconfigs[0].getAvailableRAM();
        uint64_t krambytes = divideRounded(min(args.K, vcfdata.getNReferenceHapsMax() + (args.iters > 1 ? vcfdata.getNTarget()*2 : 0)), (size_t)512)*512/8; // number of 512bit PBWT RAM words for one site (either incons or refs!) in bytes
        maxpbwtsites = avmem / (2*krambytes); // number of sites that can be buffered in FPGA RAM (for fwd and bck run)
        cout << "INFO: For this run the FPGA supports in average up to " << (maxpbwtsites/2)/fpgaconfigs[0].getNumPipelines() << " split sites per target chunk." << endl;
    }

    vector<fp_type> totconfidences;
    vector<size_t> ncalls;
    if (!args.createQuickRef) {
        totconfidences.resize(vcfdata.getNTarget(), 0.0); // for storing the phasing confidences
        ncalls.resize(vcfdata.getNTarget(), 0);           // how many variants were phased for each target (for calculating the average confidence)
    }

    // set total number of steps for StatusFile
    unsigned steps;
    if (args.createQuickRef)
        steps = 1;
    else {
        steps = args.skipImputation ? 1 : 2; // reading + imputation steps
        steps += args.iters; // phasing iterations (even if phasing is skipped, we count this step)
    }
    steps *= vcfdata.getNChunks(); // all steps repeated for each chunk
    StatusFile::updateSteps(0,steps);

    // DEBUG
    MyMalloc::printSummary(string("before chunks"));
//    ofstream memofs(args.outPrefix + ".memmap_bc");
//    MyMalloc::dumpMemMap(memofs);
//    memofs.close();
    // __DEBUG

    // process data in chunks
    for (int chunk = 0; chunk < vcfdata.getNChunks(); chunk++) {

        vcfdata.processNextChunk();

        // continue to next chunk, if the task is to create the quick reference only
        // (now, for creating the Qref only one chunk is created, so, in fact, we stop here.
        //  however, in future, the Qref may also be created in chunks, why we write "continue" here instead of "break")
        if (args.createQuickRef) {
            continue;
        }

        size_t pcapacity = roundToMultiple(vcfdata.getNSNPs(), UNITWORDS*8*sizeof(BooleanVector::data_type)) / 8;

        BooleanVector::data_type *pdata = (BooleanVector::data_type*) MyMalloc::malloc(2*vcfdata.getNTarget()*pcapacity, string("pdata_c")+to_string(chunk)); // pre-initialization done below
        vector<BooleanVector> phasedTargets(2*vcfdata.getNTarget(), BooleanVector(pdata, 2*vcfdata.getNTarget()*pcapacity, 0, toBool(Haplotype::Ref)));
        {
            auto currdata = pdata;
            for (auto &t : phasedTargets) {
                t.setData(currdata, pcapacity, 0);
                currdata += pcapacity/sizeof(BooleanVector::data_type);
            }
        }
        vector<vector<float>> phasedDosages(2*vcfdata.getNTarget());

        { // phasing block
            // TODO adjust hist factor (if set to 0) according to snprate and hetrate
            if (args.hist <= 0.0 && args.hist >= 0.0) // due to "unsafe" fp-comparison...
                args.hist = 1.0;

            if (args.skipPhasing)
                cout << "\nParsing target phases:" << endl;
            else
                cout << "\nPhasing:" << endl;

            int lockfd = 0;
            if (!usefpga && !lockfile_plain.empty() && !args.skipPhasing) {
                Stopwatch swlockp("Wait lock (p)");
                stringstream ss;
                ss << "Phasing (Chunk " << chunk+1 << "/" << vcfdata.getNChunks() << "): Waiting";
                StatusFile::updateStatus(0, ss.str());
                cout << "Waiting for CPU lock... " << flush;
                lockfd = setLock(lockdir, lockfile_plain);
                cout << "got it." << endl;
                swlockp.stop();
            }

            Hybridsys hysys(fpgas, gpus);
            PBWTPhaser phaser(hysys, vcfdata, args.num_threads, args.K, args.iters, args.expectIBDcM, args.hist, args.pErr, args.bplimitSq ? args.pErr * args.pErr : args.pErr,
                    args.impMissing, args.doPrePhasing, args.noRevPhasing, args.skipPhasing, args.debug);

            Stopwatch swp("Phasing");

            if (usefpga) {
                // This commented check is wrong since the current M is NOT the number of split sites!
                // The limit has already been checked by creating the chunks with an estimated number of split sites of M/3.
//                // check if current M works for FPGA
//                if (fpgaconfigs[0].getMaxSites() < vcfdata.getNSNPs()) {
//                    StatusFile::addError("FPGA does not support M=" + to_string(vcfdata.getNSNPs())
//                            + " sites. Maximum: " + to_string(fpgaconfigs[0].getMaxSites()));
//                    exit(EXIT_FAILURE);
//                }
                // if not set explicitly numFPGAs+2 buffers are prepared for FPGA<->host communication
                unsigned buffersFPGA = args.get<unsigned>("buffers-FPGA");
                size_t buffersizeFPGA = args.get<size_t>("buffer-size-FPGA");
                phaser.setFPGAParameters(buffersFPGA == 0 ? hysys.getFPGAs().size() * 4 : buffersFPGA, buffersizeFPGA, chrono::milliseconds(args.get<unsigned long>("timeout")), maxpbwtsites);
                term_wait = true; // need proper shutdown procedure when using FPGAs
            }
            if (usegpu) {
                // if not set explicitly numGPUs+2 buffers are prepared for GPU<->host communication
                unsigned buffersGPU = args.get<unsigned>("buffers-GPU");
                size_t buffersizeGPU = args.get<size_t>("buffer-size-GPU");
                phaser.setGPUParameters(buffersGPU == 0 ? hysys.getGPUs().size() + 2 : buffersGPU, buffersizeGPU);
            }

            phaser.phase(phasedTargets, phasedDosages, totconfidences, ncalls, chunk);

            if (usefpga)
                term_wait = false; // when using FPGAs, the critical part for user termination is over here

            // if we had a user termination request, this is the place to terminate the complete program
            if (term_request)
                exit(EXIT_FAILURE);

            // determine start phases for the next chunk
            if (chunk < vcfdata.getNChunks()-1)
                vcfdata.determineStartPhases(phasedTargets);

            swp.stop();
            if (!usefpga && !lockfile_plain.empty() && !args.skipPhasing)
                releaseLock(lockfd);


            if (args.outputPhased) {
                Stopwatch swwp("Write Phased");
                cout << "Writing phased targets..." << endl;
                vcfdata.writeVCFPhased(phasedTargets);
                swwp.stop();
            }

            // write phasing confidences after last chunk
            if (chunk == vcfdata.getNChunks()-1 && !args.skipPhasing) {
                vcfdata.writePhasedConfidences(totconfidences, ncalls);
            }

        } // END phasing block

        if (!args.skipImputation) { // imputation block

            // do reference imputation
            cout << "\nImputation:" << endl;

            stringstream ss;
            ss << "Imputation (Chunk " << chunk+1 << "/" << vcfdata.getNChunks() << ")";

            int lockfd = 0;
            if (!lockfile_plain.empty()) {
                Stopwatch swlockimp("Wait lock (i)");
                StatusFile::updateStatus(0, ss.str() + ": Waiting");
                cout << "Waiting for CPU lock... " << flush;
                lockfd = setLock(lockdir, lockfile_plain);
                cout << "got it." << endl;
                swlockimp.stop();
            }

            StatusFile::updateStatus(0, ss.str());

            Stopwatch swimp("Imputation+Write");

            unsigned num_workers = max(1u, args.num_threads - args.num_files);
            // get bunch size in number of variants
            size_t bunchsize = vcfdata.getBunchSize();
            // test if bunch size is appropriate with current chunk
            if (bunchsize > vcfdata.getNSNPsFullRef()) { // don't exceed the total number of reference variants in current chunk
                bunchsize = vcfdata.getNSNPsFullRef();
                // only set to multiple of num_workers if the memory increase would not be too much!
                if (bunchsize > 4*num_workers)
                    bunchsize = roundToMultiple(bunchsize, (size_t)num_workers);
            }
            size_t icapacity = roundToMultiple(bunchsize, UNITWORDS*8*sizeof(BooleanVector::data_type)) / 8;

            BooleanVector::data_type *idata = (BooleanVector::data_type*) MyMalloc::malloc(2*vcfdata.getNTarget()*icapacity, string("idata_c")+to_string(chunk)); // pre-initialization below
            vector<BooleanVector> imputedTargets(2*vcfdata.getNTarget(), BooleanVector(idata, 2*vcfdata.getNTarget()*icapacity, 0, toBool(Haplotype::Ref)));
            {
                auto currdata = idata;
                for (auto &t : imputedTargets) {
                    t.setData(currdata, icapacity, 0);
                    currdata += icapacity/sizeof(BooleanVector::data_type);
                }
            }
            vector<vector<float>> imputedDosages(2*vcfdata.getNTarget(), vector<float>(bunchsize));

            PBWTImputer imputer(vcfdata, args.statfile, args.num_threads, args.num_files, args.setmaxerrors, args.overwriteCalls || args.improveCalls, args.debug); // create with all available threads
            imputer.prepareImputation(phasedTargets);

            cout << "Using " << num_workers << " threads and " << args.num_files << " temporary files for imputation." << endl;
            imputer.setNumThreads(num_workers); // after finding the set-maximal matches, share available threads with writers
            vcfdata.writeVCFImputedPrepare(bunchsize);

            // alternating imputation and writing of bunches
            size_t nbunches = divideRounded(vcfdata.getNSNPsFullRef(), bunchsize * args.num_files);
            size_t startidx = 0;
            cout << "Imputing: 0%" << flush;
            int pgb = 0; // for progress bar
            for (size_t bunch = 0; bunch < nbunches; bunch++, startidx += bunchsize) {
                if (pgb == 3) { // print % after three dots
                    cout << (100*bunch/nbunches) << "%" << flush;
                    pgb = 0;
                } else {
                    cout << "." << flush;
                    pgb++;
                }
                //cout << "Imputing bunch " << bunch+1 << "/" << nbunches << " ..." << endl;
//                float progress = (bunch/(float)nbunches)/(float)vcfdata.getNChunks();
//                progress += chunk/(float)vcfdata.getNChunks();
                StatusFile::updateStatus(bunch/(float)nbunches);

                for (unsigned f = 0; f < args.num_files; f++) {
                    Stopwatch swimpb("imputeBunch");
                    imputer.imputeBunch(f, bunchsize, imputedTargets, imputedDosages);
                    swimpb.stop();
                    vcfdata.writeVCFImputedBunch(f, startidx, bunchsize, imputedTargets, imputedDosages, phasedTargets, phasedDosages);
                    // cleanup data for next block
                    memset(idata, 0, 2*vcfdata.getNTarget()*icapacity);
                }
            }
            vcfdata.writeVCFImputedClose(); // temporary files are written, data is cleaned up

            MyMalloc::free(idata);

            if (pgb == 0) // just printed "xx%"
                cout << ".";
            cout << "100%" << endl;

            // release CPU lock
            if (!lockfile_plain.empty())
                releaseLock(lockfd);

            // concat temporary files
            vcfdata.writeVCFImputedConcat();

            swimp.stop();

            StatusFile::nextStep();

        } // END imputation block

        MyMalloc::free(pdata);

        // DEBUG
        MyMalloc::printSummary(string("intermediate after chunk ")+to_string(chunk+1));
//        ofstream ofs(args.outPrefix + ".memmap_c" + to_string(chunk));
//        MyMalloc::dumpMemMap(ofs);
//        ofs.close();
        // __DEBUG

    } // end for all chunks
    StatusFile::clearContext();

    // combine output files
    if (!args.skipImputation)
        vcfdata.combineChunks();

    // print a summary of the analysis
    if (!args.createQuickRef)
        vcfdata.printSummary();

    StatusFile::updateStatus(0, "Finished");

    cout << endl;

    // dump runtime summary
    Stopwatch::dump(cerr);

    } // all destructors called

    // DEBUG
    MyMalloc::printSummary(string("at end"));
//    ofstream ofs(args.outPrefix + ".memmap_end");
//    MyMalloc::dumpMemMap(ofs);
//    ofs.close();
    // __DEBUG

    return EXIT_SUCCESS;
}
