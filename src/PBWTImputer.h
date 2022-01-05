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

#ifndef PBWTIMPUTER_H_
#define PBWTIMPUTER_H_

#include "Datatypes.h"
#include "VCFData.h"
#include "PBWT.h"
#include "TargetImp.h"

using namespace std;

class PBWTImputer {

public:
    PBWTImputer(const VCFData& vcfdata, const string &statfile, unsigned num_threads, unsigned num_blocks, int setmaxerrors, bool imputeCalls, bool debug);
    ~PBWTImputer() { if(pbwt) delete pbwt; }

    void setNumThreads(unsigned num_threads_) { num_threads = num_threads_; }
    void prepareImputation(const vector<BooleanVector> &phasedTargets);

    void imputeBunch(unsigned block, size_t nbunch, vector<BooleanVector> &imputedTargets, vector<vector<float>> &imputedDosages);

private:

    const VCFData &vcfdata;

    size_t nTargetHap; // total number of target haplotypes
    size_t K; // number of reference haplotypes for PBWT
    size_t M; // number of common sites in reference and targets

    vector<TargetImp> targets;

    PBWT *pbwt = 0;

    const string &statfile;

    unsigned num_threads;
    unsigned num_blocks;

    int setmaxerrors;
    bool imputeCalls; // explicitly do call site imputation, makes only sense with setmaxerrors > 0

    bool debug;

};

#endif /* PBWTIMPUTER_H_ */
