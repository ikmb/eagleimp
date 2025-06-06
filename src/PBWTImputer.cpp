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

#include <list>
#include <omp.h>

#include "PBWTInterval.h"
#include "Stopwatch.h"

#include "MyMalloc.h"

#include "PBWTImputer.h"

PBWTImputer::PBWTImputer(const VCFData& vcfdata_, const string &statfile_, unsigned num_threads_, unsigned num_blocks_, int setmaxerrors_, bool imputeCalls_, bool debug_) :
    vcfdata(vcfdata_),
    nTargetHap(vcfdata_.getNTarget()*2-vcfdata_.getNHaploidsTgt()),
    K(vcfdata_.getNReferenceHaps()-vcfdata_.getNHaploidsRef()),
    M(vcfdata_.getNSNPs()),
    statfile(statfile_),
    num_threads(num_threads_),
    num_blocks(num_blocks_),
    setmaxerrors(setmaxerrors_),
    imputeCalls(imputeCalls_),
    debug(debug_) {

    // generate warning if the missingness rate for some targets is too high (>50%)
    vector<size_t> missgt50tgts;
    for (size_t t = 0; t < vcfdata_.getNTarget(); t++) {
        if (vcfdata.getTgtMissPos()[t].size() > M/2)
            missgt50tgts.push_back(t);
    }
    if (missgt50tgts.size()) {
        cout << endl;
        stringstream ss;
        ss << "Missingness rate is >50% for the following targets (idx / rate / ID). Imputation will generate dosage conflicts here.\n";
        for (size_t t : missgt50tgts) {
            float mr = vcfdata.getTgtMissPos()[t].size() / (float) M;
            ss << "  " << t << " / ";
            ss << fixed << setprecision(1) << roundf(1000*mr)/10.0 << "% / ";
            ss << vcfdata.getTargetIDs()[t] << "\n";
        }
        StatusFile::addWarning(ss.str());
    }
}

void PBWTImputer::prepareImputation(const vector<BooleanVector> &phasedTargets, const vector<size_t>& num_sites_per_block) {

    // create transposed reference (reduced to target sites):
    // (only required when haploids are present)
    if (vcfdata.getNHaploidsRef()) {

        // allocate memory and create vector structure for transposed reference
        size_t capacity = roundToMultiple(K, UNITWORDS * sizeof(BooleanVector::data_type) * 8) / 8;
        BooleanVector::data_type *refdata = (BooleanVector::data_type*) MyMalloc::calloc(M * capacity, 1, string("refdata")); // initialized with 'false'

        // transposed reference reduced to target sites (used as base for PBWT)
//        vector<BooleanVector> refT = vector<BooleanVector>(M, BooleanVector(refdata, M * capacity, 0, false)); // init complete area with false
        RingBuffer<BooleanVector> refT(M,M);

        auto curr_data = refdata;
        for (auto &ref : refT) {
            ref.setData(curr_data, capacity, K);
            curr_data += capacity / sizeof(BooleanVector::data_type); // jump to next ref
        }

        // copy transpose reference data
        Stopwatch swrefT("refT");
        size_t Nref = vcfdata.getNReferenceHaps()/2;
        const auto &crefT = vcfdata.getReferenceT_const();

        omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(static, 1)
        for (size_t m = 0; m < M; m++) { // vertical
            auto &rT = refT[m];
            for (size_t n = 0, nr = 0; nr < Nref; nr++) { // horizontal write
                if (vcfdata.getHaploidsRef()[nr]) { // haploid
                    rT.setWithPreInit(n, crefT[m][2*nr]);
                    n++;
                } else {
//                    ref.setWithPreInit(n, cref[2*nr][m]);
//                    ref.setWithPreInit(n+1, cref[2*nr+1][m]);
                    rT.setPairWithPreInit(n, crefT[m][2*nr], crefT[m][2*nr+1]);
                    n+=2;
                }
            }
        }
        swrefT.stop();

        // build PBWT of the transposed reference:
        // since we want to parallelize over several targets all using the same PBWT,
        // we need to finish the complete PBWT and use it as a constant in the subroutines.
        Stopwatch swrefpbwt("refPBWT");
        pbwt = new PBWT(&refT, true); // also storing absolute permutation arrays
        pbwt->advanceTo(M-1);
        swrefpbwt.stop();

        MyMalloc::free(refdata); // not required anymore

    } else { // no haploids: directly take transposed reference that was created when reading the data
        // build PBWT of the transposed reference:
        Stopwatch swrefpbwt("refPBWT");
        pbwt = new PBWT(&(vcfdata.getReferenceT_const()), true); // also storing absolute permutation arrays
        pbwt->advanceTo(M-1);
        swrefpbwt.stop();
    }


    // Find set-maximal matches:

    Stopwatch swfindsetmax("findSetMax");
    targets.reserve(nTargetHap);
    for (size_t nhap = 0; nhap < nTargetHap; nhap++)
        targets.emplace_back(nhap, vcfdata, *pbwt, phasedTargets[vcfdata.getHaploidsTgtMap()[nhap]], num_blocks, setmaxerrors, imputeCalls, debug);

    cout << "Using " <<
    #if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT //|| defined STOPWATCH
                     1
    #else
                    num_threads
    #endif
                    << " threads for calculating set-maximal matches." << endl;

    // parallelize match finding over all targets, using the same (constant) PBWT for the reference
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
    for (size_t nhap = 2*DEBUG_TARGET_START; nhap <= 2*DEBUG_TARGET_STOP+1; nhap++)
#else
#ifndef STOPWATCH
    omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (size_t nhap = 0; nhap < nTargetHap; nhap++)
#endif
    {
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
        cout << "\nTarget: " << nhap/2 << "." << nhap%2 << endl;
#endif
        //imputeTarget(nhap, phasedTargets[nhap], pbwt, imputedTargets[nhap], imputedDosages[nhap]);
        targets[nhap].calcSMMatches();
        targets[nhap].prepareBlocks(num_sites_per_block);
    }
    swfindsetmax.stop();

    // DEBUG
#ifdef COUNTSMMATCHES
    size_t totalmatches = 0;
    size_t totalsites = 0;
    size_t total_reducable_sites = 0;
    for (size_t nt = 0; nt < nTargetHap/2; nt++) {
        for (int i = 0; i < 2; i++) {
            cout << nt << "." << i << ":\t" <<
                    targets[2*nt + i].num_matches << " matches.\t" <<
                    targets[2*nt + i].num_sites << " sites (" << targets[2*nt+i].reducable_sites << " reducable).\t" <<
                    "avl: " << targets[2*nt + i].av_match_length << "\t" <<
                    "avh: " << targets[2*nt + i].av_match_height << endl;
            totalmatches += targets[2*nt + i].num_matches;
            totalsites += targets[2*nt + i].num_sites;
            total_reducable_sites += targets[2*nt + i].reducable_sites;
            // more DEBUG - print matches
            targets[2*nt+i].printSMMatches();
        }
    }
    cout << "Total matches: " << totalmatches << endl;
    cout << "Total sites: " << totalsites << " (" << total_reducable_sites << " removable / " << (total_reducable_sites*100/(double)totalsites) << "%)" << endl;
    cout << "M: " << M << endl;
#endif
    // __DEBUG

}

void PBWTImputer::imputeBunch(unsigned block, size_t bsize, vector<BooleanVector> &imputedTargets, vector<vector<float>> &imputedDosages) {
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
    cout << "Bunch size: " << bsize << endl;
    for (size_t nhap = 2*DEBUG_TARGET_START; nhap <= 2*DEBUG_TARGET_STOP+1; nhap++)
#else
#ifndef STOPWATCH
    omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (size_t nhap = 0; nhap < nTargetHap; nhap++)
#endif
    {
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT// || defined DEBUG_TARGET_SILENT
        cout << "Target: " << nhap/2 << "." << nhap%2 << endl;
#endif
        targets[nhap].imputeBunch(block, bsize, imputedTargets[vcfdata.getHaploidsTgtMap()[nhap]], imputedDosages[vcfdata.getHaploidsTgtMap()[nhap]]);
    }
}


