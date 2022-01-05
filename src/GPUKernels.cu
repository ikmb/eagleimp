#ifdef USE_CUDA_GPU

#include "GPUEngine.h"
#include "GPUKernels.h"

// global data information
__constant__ __device__ unsigned devK = 0;
__constant__ __device__ unsigned devKwords = 0;
__constant__ __device__ unsigned devM = 0;
__constant__ __device__ unsigned devNrefhapsAD = 0;
__constant__ __device__ size_t devNrefhapsADwords = 0;
__constant__ __device__ size_t devTgtGTSize = 0;
__constant__ __device__ size_t devTgtSplitSize = 0;
__constant__ __device__ size_t devTgtBHSize = 0;

//// shared memory
//extern __shared__ unsigned char sharedmem[];

// host wrapper for copying the constant symbols
__host__ void copyConstantsToDevice(
        unsigned K,
        unsigned Kwords,
        unsigned M,
        unsigned NrefhapsAD,
        size_t NrefhapsADwords,
        size_t tgtGTSize,
        size_t tgtSplitSize,
        size_t tgtBHSize) {
    checkCUDAError(cudaMemcpyToSymbol(devK, &K, sizeof(unsigned), 0, cudaMemcpyHostToDevice))
    checkCUDAError(cudaMemcpyToSymbol(devKwords, &Kwords, sizeof(unsigned), 0, cudaMemcpyHostToDevice))
    checkCUDAError(cudaMemcpyToSymbol(devM, &M, sizeof(unsigned), 0, cudaMemcpyHostToDevice))
    checkCUDAError(cudaMemcpyToSymbol(devNrefhapsAD, &NrefhapsAD, sizeof(unsigned), 0, cudaMemcpyHostToDevice))
    checkCUDAError(cudaMemcpyToSymbol(devNrefhapsADwords, &NrefhapsADwords, sizeof(size_t), 0, cudaMemcpyHostToDevice))
    checkCUDAError(cudaMemcpyToSymbol(devTgtGTSize, &tgtGTSize, sizeof(size_t), 0, cudaMemcpyHostToDevice))
    checkCUDAError(cudaMemcpyToSymbol(devTgtSplitSize, &tgtSplitSize, sizeof(size_t), 0, cudaMemcpyHostToDevice))
    checkCUDAError(cudaMemcpyToSymbol(devTgtBHSize, &tgtBHSize, sizeof(size_t), 0, cudaMemcpyHostToDevice))
}

//__device__ int countZeros(const crefword_type* data) {
//    int cnt1 = 0;
//    if (sizeof(crefword_type) == 8) {
//        for (uint32_t i = 0; i < devKwords; i++) {
//            cnt1 += __popcll(((const unsigned long long*)data)[i]); // CUDA popcnt intrinsic (64bit)
//        }
//    } else { // assuming crefword_type is <= 4 here...
//        for (uint32_t i = 0; i < devKwords; i++) {
//            cnt1 += __popc((unsigned) (data[i])); // CUDA popcnt intrinsic (32bit)
//        }
//    }
//    // padding to full words is assumed to be 0 and not counted!!! Thus, cnt0 = K - cnt1;
//    return devK-cnt1;
//}

__device__ inline void mergeInconsistencies(data_type* dest, const data_type* src, bool gtis0, bool gtis2, unsigned words) {
//    unsigned thrptgt = blockDim.x;
//    unsigned thr = threadIdx.x;
    // determine inconsistency according to the current genotype
    if (gtis0 && !gtis2) { // HomRef
        // all 1-haps are inconsistent to the 0-genotype
        for (unsigned i = 0; i < words; i++) {
            dest[i] |= src[i];
        }
    } else if (!gtis0 && gtis2) { // HomAlt
        // all 0-haps are inconsistent to the 1-genotype
        for (unsigned i = 0; i < words; i++) {
            dest[i] |= ~(src[i]);
        }
    }
    // else: missing or heterozygous genotypes are consistent to everything, so do nothing
}

__device__ inline void copyAll(data_type* dest, const data_type* src, int* count0) {
    memcpy(dest, src, devNrefhapsADwords * sizeof(data_type)); // need to copy including padding since destination was not cleared before
    int count1 = 0;
    for (unsigned i = 0; i < (unsigned)devNrefhapsADwords; i++) {
        count1 += __popcll(src[i]);
    }
    *count0 = devK - count1;
}

__device__ inline void copyBestHaps(data_type* dest, const data_type* src, const data_type* besthaps, int* count0) {
    // it is assumed that the destination area was cleared before
    data_type* curr_wrptr = dest;
    data_type curr_wr = 0ull;
    data_type curr_wrmask = 1ull;
    const unsigned* curr_rdptr = (const unsigned*)src;
    unsigned curr_rd = *curr_rdptr;
    const unsigned* curr_bhptr = (const unsigned*)besthaps;
    unsigned curr_bh = *curr_bhptr;

    unsigned currhap = 0;
    uint32_t count1 = 0;
    while (currhap < devNrefhapsAD) {
        if (curr_bh) { // not zero: find first set bit (ffs)
            int jump = __ffs(curr_bh);
            curr_rd >>= jump-1; // jump over not required bits (zero-bits), LSB is now the current bit to read
            if (curr_rd & 1) { // bit is 1 -> copy!
                curr_wr |= curr_wrmask;
                count1++;
            }
            curr_wrmask <<= 1; // shift to next write position
            if (curr_wrmask == 0ull) {
                curr_wrmask = 1ull;
                *curr_wrptr = curr_wr;
                curr_wr = 0ull;
                curr_wrptr++;
            }
            curr_rd >>= 1; // next potential read position
            if (jump == 32) // to avoid undefined behavior
                curr_bh = 0;
            else
                curr_bh >>= jump; // next potential best hap
        } else { // no best haps in this word (anymore): take next word
            curr_rdptr++;
            curr_rd = *curr_rdptr;
            curr_bhptr++;
            curr_bh = *curr_bhptr;
            currhap += 32;
        }
    }
    if (curr_wrmask > 1ull)
        *curr_wrptr = curr_wr;

    *count0 = devK - count1;
}

__global__ void crefTKernel(
        const data_type* devRefTData,
        const uint32_t* devCrefTOffsets,
        const data_type* devTgtData,
        const data_type* devBHData,
        data_type* devCrefTData,
        data_type* devIncFullData,
        int* devCount0Data,
        int ntargets) {

    uint32_t gid = threadIdx.x + blockIdx.x * blockDim.x; /*global id*/
//    if (gid >= ntargets)
//        return;

    unsigned tgt = blockIdx.x;
    unsigned thrptgt = blockDim.x;
    unsigned thr = threadIdx.x;
//    unsigned thrptgt = 16;
//    unsigned tgt = blockIdx.x/thrptgt;
//    unsigned thr = blockIdx.x%thrptgt;

    if (tgt >= ntargets)
        return;

//    printf("tgt: %d thr: %d (thrptgt: %d)\n", tgt, thr, thrptgt);

    // source
    const data_type* curr_tgtgtsptr = devTgtData + tgt*((devTgtGTSize+devTgtSplitSize)/sizeof(data_type));
    data_type curr_tgtgts0 = *curr_tgtgtsptr;
    data_type curr_tgtgts2 = *(curr_tgtgtsptr+1);
    const data_type* curr_splitsptr = curr_tgtgtsptr + devTgtGTSize/sizeof(data_type);
    data_type curr_splits = *curr_splitsptr;
    data_type curr_mask = 1ull;

    const data_type* besthaps = devBHData + tgt*(devTgtBHSize/sizeof(data_type));

    // reference
    const data_type* curr_refhaps = devRefTData;

    // destination
    data_type* incfulldata = devIncFullData + gid * devNrefhapsADwords;
    memset(incfulldata, 0, devNrefhapsADwords * sizeof(data_type));

    data_type* currrefincTptr = devCrefTData + devCrefTOffsets[tgt]*devKwords;
    if (thr == 0 && devK != devNrefhapsAD) // clear destination area only if the data is not going to be copied as whole anyway
        memset(currrefincTptr, 0, (devCrefTOffsets[tgt+1]-devCrefTOffsets[tgt])*devKwords*sizeof(data_type));

    int* currcount0ptr = devCount0Data + devCrefTOffsets[tgt];

    __syncthreads();

//    printf("id: %d tgt: %016llx spl: %016llx bh: %016llx ref: %016llx cref: %016llx if: %016llx\n", gid, (size_t)curr_tgtgts, (size_t)curr_splits, (size_t)besthaps, (size_t)curr_refhaps, (size_t) currrefincTptr, incfulldata);
//    if (gid == 0) {
//        printf("id: %d devK: %u devKwords: %u devM: %u devNrefhapsAD: %u devNrefhapsADwords: %u\n", gid, devK, devKwords, devM, devNrefhapsAD, devNrefhapsADwords);
//        for (int nt = 0; nt <= ntargets; nt++)
//            printf("off%d: %llu\n", nt, devCrefTOffsets[nt]);
//    }

    unsigned nrefhapsADwordsnopad = (devNrefhapsAD+8*sizeof(data_type)-1)/(8*sizeof(data_type));
    data_type nrefhapsADpadmask = (devNrefhapsAD % (8*sizeof(data_type))) ? ((1ull << (devNrefhapsAD % (8*sizeof(data_type)))) - 1) : ~(0ull);

    // TODO could better distribute the load by distributing over M, but need a mapping for s then.
    unsigned s = 0;
    for (unsigned m = 0; m < devM; m++) { // for all variants

        if (curr_splits & curr_mask) { // reached split site
            if (s % thrptgt == thr) {
                // copy current incon data
                if (devK == devNrefhapsAD) { // copy all data directly
                    // need to mask potentially false set padding bits
                    incfulldata[nrefhapsADwordsnopad-1] &= nrefhapsADpadmask;
                    copyAll(currrefincTptr, incfulldata, currcount0ptr);
                } else { // K is a subset
                    copyBestHaps(currrefincTptr, incfulldata, besthaps, currcount0ptr);
//                    if (tgt == 11 && s == 3) {
//                        printf("%016llx %016llx %016llx %016llx\n", *currrefincTptr, *(currrefincTptr+1), *(currrefincTptr+2), *(currrefincTptr+3));
////                        printf("%016llx %016llx %016llx %016llx\n", *incfulldata, *(incfulldata+1), *(incfulldata+2), *(incfulldata+3));
//                    }
                }

                currrefincTptr += devKwords; // points to ref (odd) part now
                currcount0ptr++;

                // clean incon data for next segment
                memset(incfulldata, 0, devNrefhapsADwords * sizeof(data_type));

                // copy current reference data
                if (devK == devNrefhapsAD) { // copy all data directly
                    copyAll(currrefincTptr, curr_refhaps, currcount0ptr);
                } else { // K is a subset
                    copyBestHaps(currrefincTptr, curr_refhaps, besthaps, currcount0ptr);
                }

                currrefincTptr += devKwords; // points to incon (even) part now
                currcount0ptr++;
            } else {
                currrefincTptr += 2*devKwords;
                currcount0ptr += 2;
            }
            s++; // count splits
        } else { // site between splits
            if (s % thrptgt == thr) {
                // generate new incon data and combine with current
                mergeInconsistencies(incfulldata, curr_refhaps, curr_tgtgts0 & curr_mask, curr_tgtgts2 & curr_mask, nrefhapsADwordsnopad);
            }
        }

        // increment all pointers
        curr_mask <<= 1;
        if (curr_mask == 0) { // means the 1 was shifted out
            curr_mask = 1ull;
            curr_tgtgtsptr += 2; // need to jump two words since is0 and is2 are stored alternately
            curr_tgtgts0 = *curr_tgtgtsptr;
            curr_tgtgts2 = *(curr_tgtgtsptr+1);
            curr_splitsptr++;
            curr_splits = *curr_splitsptr;
        }
        curr_refhaps += devNrefhapsADwords;

    } // END for all variants

    if (s % thrptgt == thr) {
        // copy incon data from last segment
        if (devK == devNrefhapsAD) { // copy all data directly
            // need to mask potentially false set padding bits
            incfulldata[nrefhapsADwordsnopad-1] &= nrefhapsADpadmask;
            copyAll(currrefincTptr, incfulldata, currcount0ptr);
        } else { // K is a subset
            copyBestHaps(currrefincTptr, incfulldata, besthaps, currcount0ptr);
        }
    }

}

// parallelization over K (DEPRECATED!)
__global__ void pbwtKernel(
        const uint32_t* devCrefTOffsets,
        const data_type* devCrefTData,
        const int* devCount0Data,
        uint32_t* devPBWTData,
        uint32_t* devPermData,
        int ntargets,
        bool reverse) {

    unsigned tgt = blockIdx.x;
    unsigned thrptgt = blockDim.x;
    unsigned thr = threadIdx.x;

    if (tgt >= ntargets)
        return;

    // number of PBWT sites for this target
    int Mpbwt = devCrefTOffsets[tgt+1] - devCrefTOffsets[tgt];

    // determine chunksize for threads, equally distributed and overlap taken into account
    const int overlap = 96;
    int totalsize = Mpbwt + (thrptgt-1)*overlap; // total number of sites including overlap
    int chunksize = totalsize/thrptgt + (totalsize%thrptgt > thr ? 1 : 0); // number of sites for each thread
    int mstart = thr * (totalsize/thrptgt) + min(totalsize%thrptgt, thr); // start offset as if all sites including overlaps do not overlap
    mstart -= thr * overlap; // remove all preliminary overlaps to get real start site

    // data pointers
    uint32_t *compr_ptr = devPBWTData + devCrefTOffsets[tgt]*devKwords*4;
    const data_type *curr_cref_ptr = devCrefTData + devCrefTOffsets[tgt]*devKwords;
    const int* curr_cnt0_ptr = devCount0Data + devCrefTOffsets[tgt];
    uint32_t *perm = devPermData + (tgt * thrptgt + thr) * 2 * devK;
    uint32_t *nperm = perm + devK;

    // correction of data pointers due to 'reverse' and/or chunk threads
    compr_ptr += mstart * devKwords*4;
    if (reverse) {
        curr_cref_ptr += (Mpbwt-1-mstart) * devKwords;
        curr_cnt0_ptr += Mpbwt-1-mstart;
    } else {
        curr_cref_ptr += mstart * devKwords;
        curr_cnt0_ptr += mstart;
    }

//    printf("tgt: %d thr: %d (thrptgt: %d) K: %d Kwords: %d Kpad: %d Off: %d Mpbwt: %d\n", tgt, thr, thrptgt, devK, devKwords, Kpad, devCrefTOffsets[tgt], Mpbwt);
//    printf("tgt: %d thr: %d K: %d Kwords: %d Kpad: %d -- Mpbwt: %d total: %d chunk: %d mstart: %d msync: %d mstop: %d\n", tgt, thr, devK, devKwords, Kpad, Mpbwt, totalsize, chunksize, mstart, mstart + chunksize - overlap, mstart + chunksize);

    // the first absPerm is the identity
    for (uint32_t i = 0; i < devK; i++)
        perm[i] = i;

    for (int midx = mstart; midx < mstart + chunksize; midx++) {

        // wait for other threads if reached the begin of the overlap
        // (doesn't matter for the last thread, would be simply the begin of the last 'overlap' sites for it)
        if (midx == mstart + chunksize - overlap)
            __syncthreads();

        uint32_t *curr_compr_ptr = compr_ptr;
        int cnt0 = *curr_cnt0_ptr;

        int p0next = 0;
        int p1next = cnt0;
        uint32_t curr_permdata = 0;
        uint32_t curr_mask = 1;
        uint32_t *permptr = perm;
        for (int i = 0; i < devK; i++) {
            uint32_t permval = *permptr++;
            data_type crefword = *(curr_cref_ptr + permval/64);
            data_type crefmask = 1ull << (permval%64);
            if (crefword & crefmask) { // curr sample is 1

                curr_permdata |= curr_mask; // store 1 in current perm vector

                // insert index of this 1-sequence in 1 part of vector
                *(nperm+p1next) = permval;
                p1next++;

            } else { // curr sample is 0

                // insert index of this 0-sequence in 0 part of vector
                *(nperm+p0next) = permval;
                p0next++;

            }
            curr_mask <<= 1; // select next bit
            if (curr_mask == 0) { // the '1' has been shifted out of the vector
//                cerr << dec << midx << " " << i << " " << (midx * Kdiv + (i/32))  << " " << hex << setw(8) << setfill('0') << curr_permdata << endl;
                *curr_compr_ptr = curr_permdata; // store permuted vector
                curr_compr_ptr++; // now points to offset field
                *curr_compr_ptr = p0next; // store current cnt0 (which is p0next) as offset
                curr_compr_ptr++; // now points to next data field
                curr_permdata = 0; // reset permutation data
                curr_mask = 1; // reset mask
            }
        }
        if (devK%32) { // need to write the last incomplete word
//            cerr << dec << midx << " " << K-1 << " " << (midx * Kdiv + ((K-1)/32)) << " " << hex << setw(8) << setfill('0') << curr_permdata << endl;
            curr_mask = ~(curr_mask-1); // set all remaining bits to one in order to insert 1-padding
            *curr_compr_ptr = curr_permdata | curr_mask; // insert 1-padding (necessary for interval mapping!)
            curr_compr_ptr++;
            *curr_compr_ptr = p0next; // same as cnt0 now
            curr_compr_ptr++;
        }

        // next site
        compr_ptr += devKwords*4; // NOTE: devKwords refers to data_type words used for K bits. Here PBWT needs 2x the size and x2 due to uint32_t words

        // The just calculated permutation array will be the input permutation of the next step. So, we simply switch the pointers.
        uint32_t* tmp = perm;
        perm = nperm;
        nperm = tmp;

        // update pointers
        if (reverse) {
            curr_cref_ptr -= devKwords;
            curr_cnt0_ptr--;
        } else {
            curr_cref_ptr += devKwords;
            curr_cnt0_ptr++;
        }

    }

}


__global__ void pbwtPartedKernel(
        const uint32_t* devCrefTOffsets,
        const data_type* devCrefTData,
        uint32_t* devPBWTData,
        uint32_t* devPermData,
        int ntargets,
        bool reverse) {

    unsigned tgt = blockIdx.x;
    unsigned thrptgt = blockDim.x;
    unsigned thr = threadIdx.x;

    if (tgt >= ntargets)
        return;

    // number of PBWT sites for this target
    int Mpbwt = devCrefTOffsets[tgt+1] - devCrefTOffsets[tgt];

    // destination data pointers for this target
    uint32_t *comprPBWT = devPBWTData + devCrefTOffsets[tgt]*devKwords*4; // beginning of destination space for PBWT
    uint32_t *permData = devPermData + tgt * 2*devK; // beginning of temporary space for this target

    { // Part I: PBWT creation

        // division in parts -> how large is my part?
        int kpartwordstotal = (devK+63)/64; // 64bit words in total
        int kpartwords = kpartwordstotal / thrptgt + (thr < kpartwordstotal % thrptgt ? 1 : 0); // divided up to all threads per target
        int kpartwordoffset = thr * kpartwords + (thr >= kpartwordstotal % thrptgt ? kpartwordstotal % thrptgt : 0); // offset for current thread
        int kpart = kpartwords * 64;
        if (thr == thrptgt-1 && devK % 64) // only the last thread may have an incomplete last word
            kpart -= 64 - devK % 64;

        // data pointers
        const data_type *crefTData = devCrefTData + devCrefTOffsets[tgt]*devKwords; // beginning of crefTData for this target
        const data_type *crefTDatapart = crefTData + kpartwordoffset; // beginning of crefTData for this thread
        uint32_t *comprPBWTpart = comprPBWT + kpartwordoffset*4;
        uint32_t *permDatapart = permData + 2*64*kpartwordoffset;

        // correction of data pointer if 'reverse'
        if (reverse)
            crefTDatapart += (Mpbwt-1) * devKwords;

        // DEBUG
//        printf("PBWT: tgt: %d thr: %d (thrptgt: %d) Mpbwt: %d K: %d Kwords: %d kpartwordstotal: %d kpartwordoffset: %d kpartwords: %d kpart: %d\n", tgt, thr, thrptgt, Mpbwt, devK, devKwords, kpartwordstotal, kpartwordoffset, kpartwords, kpart);
//        if (tgt == 0 && thr == 0)
//            printf("PBWT...\n");

        // run partial PBWT creation
        buildComprPBWTPart(crefTDatapart, comprPBWTpart, permDatapart, Mpbwt, kpart, reverse);

    } // END Part I

    __syncthreads();

    { // Part II: Merge and Fill Offsets

        // calculate pipeline depth
        unsigned totaldepth = 0;
        unsigned mydepth = 0;
        unsigned x = thrptgt / 2; // threads in current depth
        unsigned y = x; // thread index exceeding current depth
        while (x > 0) {
            totaldepth++;
            if (thr >= y)
                mydepth = totaldepth;
            x /= 2;
            y += x;
        }

        // calculate my merging part (in fact, this is also calculated by the fill thread, but not required)
        unsigned thrpdpth = 1u << (totaldepth - mydepth -1); // 2^... = number of threads that are at the same depth in the pipeline as me
        unsigned relthr = thr % thrpdpth; // relative thread index of my thread to all other threads in the same depth
        // indices of original PBWT parts for merging
        int pstart = relthr << (mydepth+1); // start index of first part
        int pmid   = (2*relthr+1) << mydepth; // start index of second part
        int pstop  = (relthr+1) << (mydepth+1); // end index of second part (exclusive)
        // as for the PBWT calculation, we calculate the offset and the sizes of each part
        int kpartwordstotal = (devK+63)/64; // the partition of the PBWT was based on 64bit words because of the original data word size
        int kpartwords = kpartwordstotal/thrptgt; // minimum number of words each thread has processed
        int koffwords = pstart * kpartwords; // minimum offset in words to start of first merging part
        koffwords += (pstart >= kpartwordstotal % thrptgt ? (kpartwordstotal % thrptgt) : (pstart % thrptgt)); // offset correction
        int k0words = pmid * kpartwords; // minimum offset in words to start of second merging part
        k0words += (pmid >= kpartwordstotal % thrptgt ? (kpartwordstotal % thrptgt) : (pmid % thrptgt)); // offset correction
        k0words -= koffwords; // actual number of words in k0 part
        int k1words = pstop * kpartwords; // minimum offset in words to (exclusive) end of second merging part
        k1words += (pstop >= kpartwordstotal % thrptgt ? (kpartwordstotal % thrptgt) : (pstop % thrptgt)); // offset correction
        k1words -= koffwords + k0words; // actual number of words in k0 part
        int koff = koffwords*64;
        int k0 = k0words*64;
        int k1 = k1words*64;
        // k's are a multiple of 64 now, but if the PBWT consists of an odd number of 32bit words (not counting the offset fields),
        // the last 32 references do not carry a correct padding, so we have to remove them here.
        if (koff+k0+k1-32 > devK)
            k1 -= 32;
        // correct temporary space by the current depth
        uint32_t *orgspace = permData + devKwords*4 * mydepth; // devKwords refers to 64bit words and we require 2x devKwords space for origins in each pipeline step

//        for (int d = 0; d < totaldepth; d++) {
//            if (d == mydepth) {
        if (thr < thrptgt-1) {
//                // DEBUG
//                if (tgt == 0 && relthr == 0)
//                    printf("d=%u merge...\n", mydepth);

                merge(comprPBWT, orgspace, Mpbwt, koff, k0, k1, mydepth, totaldepth);
        }
////            else
////                printf("thr: %u Waiting for sync.\n", threadIdx.x);
            __syncthreads();
//        }
        if (thr == thrptgt-1) { // only last thread fills the offsets
//            // DEBUG
//            if (tgt == 0)
//                printf("d=%u fill...\n", mydepth);

            fillOffsets(comprPBWT, Mpbwt, mydepth);

//            // DEBUG
//            if (tgt == 0)
//                printf("finished.\n");
        }

    } // END Part II

}

__device__ inline void buildComprPBWTPart(
        const data_type* crefTData,
        uint32_t* comprPBWT,
        uint32_t* permData,
        int Mpbwt,
        int kpart,
        bool reverse) {

    const data_type *cref_ptr = crefTData;
    uint32_t *compr_ptr = comprPBWT;
    uint32_t *perm = permData;
    uint32_t *nperm = permData + kpart;

    // the first absPerm is the identity
    for (uint32_t i = 0; i < kpart; i++)
        perm[i] = i;

    for (int midx = 0; midx < Mpbwt; midx++) {

        uint32_t *curr_compr_ptr = compr_ptr;

        int cnt1 = 0;
        const data_type *dptr = cref_ptr;
        for (int i = 0; i < (kpart+63)/64; i++)
            cnt1 += __popcll(*dptr++);

//        int p0next = 0;
//        int p1next = kpart - cnt1;
        uint32_t *nperm0 = nperm;
        uint32_t *nperm1 = nperm + kpart - cnt1;

//        // DEBUG
//        int cnt0 = p1next;
//        printf("m: %u cnt0: %u cnt1: %u\n", midx, cnt0, cnt1);

        uint32_t curr_permdata = 0;
        uint32_t curr_mask = 1;
        uint32_t *permptr = perm;
        for (int i = 0; i < kpart; i++) {
            uint32_t permval = *permptr++;
            data_type crefword = *(cref_ptr + permval/64);
            data_type crefmask = 1ull << permval%64;
            if (crefword & crefmask) { // curr sample is 1

                curr_permdata |= curr_mask; // store 1 in current perm vector

                // insert index of this 1-sequence in 1 part of vector
                *nperm1++ = permval;
//                *(nperm+p1next) = permval;
//                p1next++;

//                // DEBUG
//                if (p1next > kpart)
//                    printf("!!!");
//                printf("k: %i perm: %u wd: %016llx mask: %016llx bit: 1 p0next: %u p1next: %u\n", i, permval, crefword, crefmask, p0next, p1next);

            } else { // curr sample is 0

                // insert index of this 0-sequence in 0 part of vector
                *nperm0++ = permval;
//                *(nperm+p0next) = permval;
//                p0next++;

//                // DEBUG
//                if (p0next > cnt0)
//                    printf("!!!");
//                printf("k: %i perm: %u wd: %016llx mask: %016llx bit: 0 p0next: %u p1next: %u\n", i, permval, crefword, crefmask, p0next, p1next);

            }
            curr_mask <<= 1; // select next bit
            if (curr_mask == 0) { // the '1' has been shifted out of the vector
//                cerr << dec << midx << " " << i << " " << (midx * Kdiv + (i/32))  << " " << hex << setw(8) << setfill('0') << curr_permdata << endl;
                *curr_compr_ptr = curr_permdata; // store permuted vector
                curr_compr_ptr+=2; // jump over offset field
                curr_permdata = 0; // reset permutation data
                curr_mask = 1; // reset mask
            }
        }
        if (kpart%32) { // need to write the last incomplete word
//            cerr << dec << midx << " " << K-1 << " " << (midx * Kdiv + ((K-1)/32)) << " " << hex << setw(8) << setfill('0') << curr_permdata << endl;
            curr_mask = ~(curr_mask-1); // set all remaining bits to one in order to insert 1-padding
            *curr_compr_ptr = curr_permdata | curr_mask; // insert 1-padding (necessary for interval mapping!)
        }

        // the current count0 value is stored in the first offset field to speedup things
        compr_ptr[1] = kpart - cnt1;

        // next site
        compr_ptr += devKwords*4; // NOTE: devKwords refers to data_type words used for K bits. Here PBWT needs 2x the size and x2 due to uint32_t words

        // The just calculated permutation array will be the input permutation of the next step. So, we simply switch the pointers.
        uint32_t* tmp = perm;
        perm = nperm;
        nperm = tmp;

        // update pointers
        if (reverse) {
            cref_ptr -= devKwords;
        } else {
            cref_ptr += devKwords;
        }

    }

}

// ptr is incremented by 2 if mask gets an overflow
__device__ inline void rdincrement(uint32_t* &ptr, uint32_t &mask, uint32_t &ptrval) {
    mask <<= 1;
    if (mask == 0u) {
        mask = 1u;
        ptr+=2; // increment by 2 to jump over the offset field
        ptrval = *ptr;
    }
}

// ptr0 is incremented by 2 and ptr1 by 1 if mask gets an overflow
__device__ inline void rdincrement(uint32_t* &ptr0, uint32_t* &ptr1, uint32_t &mask, uint32_t &ptrval0, uint32_t &ptrval1) {
    mask <<= 1;
    if (mask == 0u) {
        mask = 1u;
        ptr0+=2; // increment by 2 to jump over the offset field
        ptr1++;
        ptrval0 = *ptr0;
        ptrval1 = *ptr1;
    }
}

// ptr0 is incremented by 2 and ptr1 by 1 if mask gets an overflow
__device__ inline void wrincrement(uint32_t* &ptr0, uint32_t* &ptr1, uint32_t &mask, uint32_t &ptrval0, uint32_t &ptrval1) {
    mask <<= 1;
    if (mask == 0u) {
        mask = 1u;
        *ptr0 = ptrval0;
        *ptr1 = ptrval1;
        ptr0+=2; // increment by 2 to jump over the offset field
        ptr1++;
        ptrval0 = 0u;
        ptrval1 = 0u;
    }
}

// parts to merge have to be located next to each other
// koff: beginning of parts to be merged (multiple of 32!)
// k0: length (in bits) of part 0 (multiple of 32!)
// k1: length (in bits) of part 1 (multiple of 32!)
//     NOTE: for the last part koff + k0 + k1 may exceed K, so be sure to have included the correct padding!
__device__ inline void merge(uint32_t* comprPBWT, uint32_t *orgspace, int Mpbwt, int koff, int k0, int k1, int pipedepth, int totaldepth) {
    // odd: the final cPBWT is stored as one 32bit word reference information and one 32bit word offset. The offset is only required in the final version,
    //      thus merging uses the memory space to store the merged information. if 'odd' is true, the information to be merged is assumed to be in the
    //      offset memory already and the output of the merge will be stored in the data memory. if 'odd' is false, the merge output will be written
    //      to the offset space.
    bool odd = pipedepth % 2;

//    // DEBUG
//    printf("merge: thr: %u Mpbwt=%u koff=%u k0off=%u k1off=%u k0=%u k1=%u k0+k1=%u depth=%u odd=%s\n", threadIdx.x, Mpbwt, koff, koff+k0, koff+k0+k1, k0, k1, k0+k1, pipedepth, odd ? "T" : "F");

    uint32_t *compr_ptr = comprPBWT + koff/16; // koff is a multiple of 32 anyway, and we need to jump over offset fields
    const int kw0 = k0/32; // k0 is a multiple of 32
    const int kw1 = k1/32; // k1 is a multiple of 32

//    uint32_t *rdptr = compr_ptr + (odd ? 1 : 0);
//    uint32_t *wrptr = compr_ptr + (odd ? 0 : 1);
//    uint32_t cnt0_0 = *wrptr; // cnt0 is located in the offset field, so the first write position
//    uint32_t cnt0_1 = *(wrptr + 2*kw0);
//    // copy data:
//    // first word and new count0 value
//    *wrptr = *rdptr;
//    *rdptr = cnt0_0 + cnt0_1;
//    for (int kw = 1; kw < kw0+kw1; kw++) { // k0 data and k1 data are required to be next to each other
//        // increment first because first word was already copied
//        rdptr += 2; // always jump over offset fields (which is where we just copied to (if odd == false))
//        wrptr += 2;
//        *wrptr = *rdptr;
//    }

    uint32_t cnt0_0 = 0;
    uint32_t cnt0_1 = 0;

    // origins
    uint32_t *origin = orgspace + koff/16; // koff is a multiple of 32 anyway, and we require two times kwords space for origins (actual and new)
    uint32_t *norigin = origin + kw0 + kw1;
    memset(origin, 0, kw0*sizeof(uint32_t)); // first k0 elements are from 0
    memset(origin+kw0, 0xff, kw1*sizeof(uint32_t)); // remaining k1 elements are from 1

    // first site does never change, so we leave it as is. we start the process from the second site
    for (int midx = 1; midx < Mpbwt+totaldepth; midx++) {
        int m = midx - pipedepth;

        if (m > 0 && m < Mpbwt) { // only continue if we have a valid site

            // if we are at the first site, we need to read and write cnt0 values
            if (m == 1) {
                cnt0_0 = *(compr_ptr + 1); // cnt0 is located in the offset field
                cnt0_1 = *(compr_ptr + 2*kw0 + 1);
                *(compr_ptr+1) = cnt0_0 + cnt0_1; // new combinded cnt0 value
            }

            // merged data from last site
            uint32_t *rdlastptr = compr_ptr + (odd || m == 1 ? 0 : 1); // where did I just write the merge results from previous site to?
            uint32_t rdlastval = *rdlastptr;
            uint32_t *rdorgptr = origin; // for reading the origin
            uint32_t rdorg = *rdorgptr;
            uint32_t rdlastmask = 1u;

            // increment to current site
            compr_ptr += 4*devKwords;

            // read source data from this site
            uint32_t *rdptr0_0 = compr_ptr + (odd ? 1 : 0);
            uint32_t rdmask0_0 = 1u;
            uint32_t rdval0_0 = *rdptr0_0;
            uint32_t *rdptr0_1 = compr_ptr + 2*kw0 + (odd ? 1 : 0);
            uint32_t rdmask0_1 = 1u;
            uint32_t rdval0_1 = *rdptr0_1;
            uint32_t *rdptr1_0 = rdptr0_0 + 2*(cnt0_0 / 32);
            uint32_t rdmask1_0 = 1u << (cnt0_0 % 32);
            uint32_t rdval1_0 = *rdptr1_0;
            uint32_t *rdptr1_1 = rdptr0_1 + 2*(cnt0_1 / 32);
            uint32_t rdmask1_1 = 1u << (cnt0_1 % 32);
            uint32_t rdval1_1 = *rdptr1_1;

            // write pointers
            uint32_t *wrptr0 = compr_ptr + (odd ? 0 : 1);
            uint32_t *wrorgptr0 = norigin; // for writing the origin
            uint32_t wrmask0 = 1u;
            uint32_t wrval0 = 0u;
            uint32_t wrorg0 = 0u;
            uint32_t *wrptr1 = wrptr0 + 2*((cnt0_0+cnt0_1)/32);
            uint32_t *wrorgptr1 = wrorgptr0 + ((cnt0_0+cnt0_1)/32);
            uint32_t wrmask1 = 1u << ((cnt0_0+cnt0_1)%32);
            uint32_t wrval1 = 0u;
            uint32_t wrorg1 = 0u;
            bool wroverlap01 = (cnt0_0+cnt0_1)%32; // need to pay attention if the last 0 word contains 1 data as well
            // read cnt0 values for this site
            uint32_t ncnt0_0 = *wrptr0;
            uint32_t ncnt0_1 = *(wrptr0 + 2*kw0);

            // directly write new combined count0 value for this site (to where I just read the first word from)
            *rdptr0_0 = ncnt0_0 + ncnt0_1;

    //        // DEBUG
    //        cerr << "m=" << midx << ": cnt0_0=" << cnt0_0 << " cnt1_0=" << (k0-cnt0_0) << " cnt0_1=" << cnt0_1 << " cnt1_1=" << (k1-cnt0_1) << endl;
    //        int wr0_0 = 0, wr0_1 = 0, wr1_0 = 0, wr1_1 = 0;

            // iterate through permutation from previous site and read accordingly from source to generate merged result for current site
            for (int k = 0; k < k0+k1; k++) {
                if (rdlastval & rdlastmask) { // last perm was 1
                    if (rdorg & rdlastmask) { // bit came from second part
                        if (rdval1_1 & rdmask1_1) { // next bit is 1
                            wrval1 |= wrmask1;
                        } else { // next bit is 0
                            ;
                        }
                        wrorg1 |= wrmask1;
                        rdincrement(rdptr1_1, rdmask1_1, rdval1_1);
    //                    wr1_1++;
                    } else { // bit came from first part
                        if (rdval1_0 & rdmask1_0) { // next bit is 1
                            wrval1 |= wrmask1;
                        } else { // next bit is 0
                            ;
                        }
                        rdincrement(rdptr1_0, rdmask1_0, rdval1_0);
    //                    wr1_0++;
                    }
                    wrincrement(wrptr1, wrorgptr1, wrmask1, wrval1, wrorg1);
                } else { // last perm was 0
                    if (rdorg & rdlastmask) { // bit came from second part
                        if (rdval0_1 & rdmask0_1) { // next bit is 1
                            wrval0 |= wrmask0;
                        } else { // next bit is 0
                            ;
                        }
                        wrorg0 |= wrmask0;
                        rdincrement(rdptr0_1, rdmask0_1, rdval0_1);
    //                    wr0_1++;
                    } else { // bit came from first part
                        if (rdval0_0 & rdmask0_0) { // next bit is 1
                            wrval0 |= wrmask0;
                        } else { // next bit is 0
                            ;
                        }
                        rdincrement(rdptr0_0, rdmask0_0, rdval0_0);
    //                    wr0_0++;
                    }
                    wrincrement(wrptr0, wrorgptr0, wrmask0, wrval0, wrorg0);
                }
                rdincrement(rdlastptr, rdorgptr, rdlastmask, rdlastval, rdorg);
            } // END for k

    //        cerr << dec << " wr0_0=" << wr0_0 << " wr1_0=" << wr1_0 << " wr0_1=" << wr0_1 << " wr1_1=" << wr1_1 << endl;

    //        // check: there should not be a reminder of 1-words or of 0-words if there is no overlap
    //        if (!wroverlap01 && wrmask0 != 1u)
    //            cerr << "!wroverlap01 && wrmask0 != 1: THIS SHOULD NOT HAPPEN!!!" << endl;
    //        if (wrmask1 != 1u)
    //            cerr << "wrmask1 != 1: THIS SHOULD NOT HAPPEN!!!" << endl;

            // write potential remainder
            if (wroverlap01) { // there is a remainder of the 0-words that has to be combined with the first 1-word
                *wrptr0 |= wrval0;
                *wrorgptr0 |= wrorg0;
            } // else: there is no remainder

            // take over cnt0 values for next site
            cnt0_0 = ncnt0_0;
            cnt0_1 = ncnt0_1;
            // switch origin vectors
            uint32_t *tmp = origin;
            origin = norigin;
            norigin = tmp;
        } // END if site is valid

        __syncthreads();

    } // END for m

//    // DEBUG
//    printf("merge: thr: %u finished.\n", threadIdx.x);
}

// generate the offset information in the merged cPBWT. if 'odd' is set, the reference information is stored in the offset field, so this function will copy
// it to the data field while calculating the offset
__device__ inline void fillOffsets(uint32_t* comprPBWT, int Mpbwt, int pipedepth) {
    bool odd = pipedepth % 2;

//    // DEBUG
//    printf("fill: thr: %u Mpbwt=%u depth=%u odd=%s\n", threadIdx.x, Mpbwt, pipedepth, odd ? "T" : "F");

    uint32_t *compr_ptr = comprPBWT;
    for (int m = 0; m < Mpbwt; m++) {
        uint32_t *curr_compr_ptr = compr_ptr;
        uint32_t cnt0 = 0;
        for (int kw = 0; kw < (devK+31)/32; kw++) {
            uint32_t data = odd && m > 0 ? *(curr_compr_ptr+1) : *curr_compr_ptr;
            cnt0 += 32-__popc(data);
            if (odd && m > 0)
                *curr_compr_ptr = data;
            curr_compr_ptr++;
            *curr_compr_ptr = cnt0;
            curr_compr_ptr++;
        }
        compr_ptr += 4*devKwords;
    }
//    // DEBUG
//    printf("fill: thr: %u finished.\n", threadIdx.x);
}

#endif // USE_CUDA_GPU


