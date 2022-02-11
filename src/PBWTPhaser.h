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

#ifndef PBWTPHASER_H_
#define PBWTPHASER_H_

#include <vector>
#include <memory>
#include <atomic>

#include "hybridsys/Hybridsys.h"

#include "Datatypes.h"
#include "VCFData.h"
#include "Target.h"

using namespace std;

class PBWTPhaser {

public:
    typedef vector<BooleanVector> refincs_type;

    typedef shared_ptr<vector<Target*>> targetQueue_type;
    struct condensedRefQueue_type {
        Target* target = NULL;
        vector<BooleanVector>* cref = NULL;
        vector<int>* count0 = NULL;
    };
    struct cPBWTQueue_type {
        Target* target = NULL;
        uint32_t* cpbwt = NULL;
        vector<int>* count0 = NULL;
    };
    struct confidence_type {
        size_t id;
        float totalconf;
        size_t ncalls;
    };

    static atomic<bool> terminate;

    PBWTPhaser(hybridsys::Hybridsys &hysys, VCFData &vcfdata, unsigned numthreads,
            size_t Karg, uint32_t iters, fp_type expectIBDcM,
            fp_type hist, fp_type pErr, fp_type pLimit, bool impMissing,
            bool doPrePhasing, bool noRevPhasing, bool skipPhasing, bool debug);

    void setFPGAParameters(unsigned num_buffers, size_t buffersize, chrono::milliseconds timeout, size_t maxrefincsites_) {
        fpga_num_buffers = num_buffers;
        fpga_buffersize = buffersize;
        fpga_timeout = timeout;
        maxpbwtsites = maxrefincsites_;
        usefpga = true;
    }

    void setGPUParameters(unsigned num_buffers, size_t buffersize) {
        gpu_num_buffers = num_buffers;
        gpu_buffersize = buffersize;
        usegpu = true;
    }

    // do the phasing!
    void phase(vector<BooleanVector> &phasedTargets, vector<vector<float>> &phasedDosages, vector<fp_type> &totconfidences, vector<size_t> &ncalls, int chunk);

private:

    // do phasing with FPGA support
    void phaseFPGA(vector<BooleanVector> &phasedTargets, vector<vector<float>> &phasedDosages, int chunk,
            vector<fp_type> &avconfidences, vector<size_t> &ncalls);

    // do phasing with GPU support
    void phaseGPU(vector<BooleanVector> &phasedTargets, vector<vector<float>> &phasedDosages, int chunk,
            vector<fp_type> &avconfidences, vector<size_t> &ncalls);

    // do phasing without additional hardware support
    void phaseCPU(vector<BooleanVector> &phasedTargets, vector<vector<float>> &phasedDosages, int chunk,
            vector<fp_type> &avconfidences, vector<size_t> &ncalls);

    unsigned numthreads;

    VCFData &vcfdata;

//    const vector<GenotypeVector> &targetsFull; // full information of targets
    size_t nTarget; // total number of target samples

    size_t Karg;
    uint32_t iters;
    fp_type expectIBDcM;
    fp_type hist; // history factor
    fp_type pErr; // switch error probability
    fp_type pLimit; // probability deviation limit in beam (usually switch error probability or squared sw err prob)

    // FPGA/GPU accelerator parameters
    hybridsys::Hybridsys &hysys;
    bool usefpga = false;
    bool usegpu = false;
    unsigned fpga_num_buffers = 0;
    size_t   fpga_buffersize = 0;
    unsigned gpu_num_buffers = 0;
    size_t   gpu_buffersize = 0;
    chrono::milliseconds fpga_timeout = chrono::milliseconds(0);
    size_t maxpbwtsites = ~0ull; // limit of PBWT sites (i.e. either ref or inc site) per target block (due to memory requirements)

    bool impMissing;
    bool doPrePhasing;
    bool doRevPhasing;
    bool skipPhasing;

    const vector<string> &targetIDs; // target ID information

    bool debug;

};

#endif /* PBWTPHASER_H_ */
