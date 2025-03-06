/*
 *    Copyright (C) 2018-2025 by Lars Wienbrandt,
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

#ifndef VCFDATA_H_
#define VCFDATA_H_

#include <tbb/concurrent_queue.h>
#include <string>
#include <vector>
#include <random>
#include <thread>
#include <fstream>

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

#include "Args.h"
#include "Datatypes.h"
#include "RingBuffer.h"
#include "MapInterpolater.h"
#include "FPGAConfigurationEagleImp.h"

#include "MyMalloc.h"

using namespace std;

static const float AF_UNKNOWN = bcf_float_missing; // used if the allele frequency is not calculated

// stupid htslib does not define this in header, so I had to copy it... :-/
#ifndef MAX_CSI_COOR
#define MAX_CSI_COOR ((1LL << (14 + 30)) - 1)
#endif

class VCFData {

public:
    VCFData(const Args &args, int argc, char**argv, const vector<FPGAConfigurationEagleImp> &fpgaconfs);

    ~VCFData() {
        if (doImputation)
            infofile.close();
        if (refdata) MyMalloc::free(refdata);
        if (refdataT) MyMalloc::free(refdataT);
        if (multiAllFlagsFullRefRegion.getData()) MyMalloc::free(multiAllFlagsFullRefRegion.getData());
        for (auto &r : referenceFullT) {
            if (r.getData()) MyMalloc::free(r.getData());
        }
        if (sr) bcf_sr_destroy(sr);
//        if (tgt_hdr) -> destroyed while destroying sr
//            bcf_hdr_destroy(tgt_hdr);
        if (tgt_hdr_cpy)
            bcf_hdr_destroy(tgt_hdr_cpy);
        if (imp_hdr)
            bcf_hdr_destroy(imp_hdr);
    }

    bool getUseFPGA() { return usefpga; }

    int getNChunks() const { return nChunks; }
    // processes the next chunk, i.e. continue reading haplotype data and prepare result data fields for phasing and imputation
    void processNextChunk();

    const vector<BooleanVector> &getReference_const() const { return reference; }
    const RingBuffer<BooleanVector> &getReferenceT_const() const { return referenceT; }
    vector<BooleanVector> &getReference() { return reference; }
    RingBuffer<BooleanVector> &getReferenceT() { return referenceT; }
    const vector<BooleanVector> &getReferenceFullT() const { return referenceFullT; }
    const vector<GenotypeVector> &getTargets() const { return targets; }
    const vector<RingBuffer<size_t>> &getTgtMissPos() const { return tgtmiss; }
    const vector<RingBufferBool> &getTgtInPhases() const { return tgtinphase; }
    const RingBuffer<fp_type> &getCMPos() const { return cMs; }
    const vector<string> &getTargetIDs() const { return targetIDs; }

    // returns the number of reference haplotypes (NOTE: haploids are encoded as homozygous diploid!)
    size_t getNReferenceHaps() const { return Nrefhaps; }
    size_t getNReferenceHapsMax() const { return Nrefhapsmax; }
    // returns the number of target samples
    size_t getNTarget() const { return Ntarget; }
    size_t getNSNPs() const { return M; }
    size_t getNSNPsFullRef() const { return Mref; }
    size_t getSNPIndexToFullRef(size_t tgtsnpindex) const { return indexToRefFull[tgtsnpindex]; }
    float getAlleleFreqFullRefRegion(size_t fullrefsnpindex) const { return alleleFreqsFullRefRegion[fullrefsnpindex]; }
    float getAlleleFreqCommon(size_t tgtsnpindex) const { return alleleFreqsCommon[tgtsnpindex]; }
    // returns a vector of flags for each sample in the reference if it is haploid (true) or not (false)
    const vector<bool> &getHaploidsRef() const { return haploidsRef; }
    const vector<bool> &getHaploidsTgt() const { return haploidsTgt; }
    const vector<size_t> &getHaploidsRefMap() const { return haploidsRefMap; }
    const vector<size_t> &getHaploidsTgtMap() const { return haploidsTgtMap; }
    // returns the number of haploid samples in the reference (without targets)
    size_t getNHaploidsRef() const { return nHaploidsRef; }
    // returns the number of haploid samples in the target
    size_t getNHaploidsTgt() const { return nHaploidsTgt; }
    unsigned getNumWorkers() const { return num_workers; }
    unsigned getNumFiles() const { return num_files; }
    // returns the number of sites per file, whereby the first and the last element indicate the number of sites in the corresponding overlap, that should be ignored for writing
    const vector<size_t>& getNumSitesPerFile() const { return num_sites_per_file; }
    size_t getBunchSize() const { return bunchsize; }
    size_t getNBunches() const { return nbunches; }

    // determines the start phases for the next chunk from the phasing results of the current chunk
    void determineStartPhases(const vector<BooleanVector> &phasedTargets);
    // return start phase for selected target, i.e. the maternal haplotype of the first heterozygous site of the overlap region of the previous chunk
    Haplotype getStartPhase(size_t targetidx) const { return toHap(chunkStartPhase[targetidx]); }

    void writePhasedConfidences(const vector<fp_type> &totalconfs, const vector<size_t> &ncalls);
    void writeVCFPhased(const vector<BooleanVector> &phasedTargets);
    void writeVCFImputedPrepare();
    void writeVCFImputedBunch(
            unsigned block,
            size_t startidx,
            size_t bsize,
            const vector<BooleanVector> &imputedTargets,
            const vector<vector<float>> &imputedDosages,
            const vector<BooleanVector> &phasedTargets,
            const vector<vector<float>> &phasedDosages);
    void writeVCFImputedClose();
    void writeVCFImputedConcat();
    void combineChunks() const;

    void printSummary() const;

    static int chrNameToNumber(const string &name, string &chrlit, bool *literally = NULL);

private:
    // reads metadata and determines the number of chunks, called in constructor
    inline void processMeta(const string &refFile, const string &vcfTarget, const string &vcfExclude, const Args &args, const vector<FPGAConfigurationEagleImp> &fpgaconfs);
    // initializes the header used for the imputation output
    inline void initImpHeader();

    // initialize variant info file
    inline void initInfoFile();
    // add excluded variant to info file
    inline void addToInfoFileExcluded(bcf1_t* tgt, const string& explanation);
    inline void addToInfoFileExcluded(size_t pos, const string& tgtID, const string& tref, const string& talt, const string& explanation);
    // add included variant to info file
    inline void addToInfoFileIncluded(size_t pos, const string& tgtID, const string& tref, const string& talt, const string& refID, const string& refall, const string& altall, const string& explanation);

    inline void processReferenceOnlySNP(bcf1_t *ref, void **ref_gt, int *mref_gt, int &numMissing, int &numUnphased, bool multiallelic);
    inline void processReferenceSNP(int nsmpl, bcf1_t *ref, void **ref_gt, int *mref_gt,
            int &numMissing, int &numUnphased, float &af, bool addForPhasing);
    inline void processTargetSNP(int nsmpl, int ngt, const int32_t *gt, bool refAltSwap, int &numMissing, size_t tgtvariantpos);
    inline void processMap(const MapInterpolater &mapInterpolater);
    inline void qRefAppendVariant(const BooleanVector& var);
    inline void qRefWriteMetaAndConcat();
    inline void qRefOpenReadMeta();
    // read all variants from current stream position until either the position is reached exactly
    // or until the closest variant BEFORE that position.
    // data is added with "push_back" to referenceFullT.
    // returns true, if at least tgt position matches, and only then, the following bools are set accordingly:
    // - if ref and alt allele were swapped, refaltswapped is set.
    // - if the reverse complements of the alleles match, strandflipped is set.
    // - if the tgt ref/alt pair has no counterpart in the reference, refalterror is set.
    // qridx is either the index (relative to region) where the tgt was found or the index of the following variant if not found.
    // NOTE: the function may swap entries at positions greater or equal to the returned qridx
    //       (happens if tgt was found in a multi-allelic variant, then the matching variant is swapped to qridx, while others are moved to higher indices)
    inline bool loadAndFindInQref(bcf1_t *tgt, size_t &qridx, bool &refaltswapped, bool &strandflipped, bool &refalterror);
    // loads the next variant from the Qref data stream and adds the data to referenceFullT by using emplace_back() if store == true
    // ( if store == false only a placeholder without data is added to referenceFullT )
    inline void qRefLoadNextVariant(bool store);

    // split multi-allelic reference SNPs:
    // *ref will be updated to a bi-allelic SNP with the reference and first alternative allele.
    // for the other alternative alleles a new record will be created with the same reference allele and put into masplits.
    // nSplitted will be incremented by the number of resulting splits (including the original).
    // gFilled are the number of haplotypes that had to be filled with random values in each split.
    inline void splitSNP(bcf1_t *ref, void **ref_gt, int *mref_gt,
            vector<bcf1_t*> &masplits, size_t &nSplitted, size_t &gFilled);

    inline void appendVersionToBCFHeader(bcf_hdr_t *hdr) const;
    inline void concatFiles(const vector<string>& filenames) const;
    inline string getOutputSuffix();
    // returns the number of variants only from the first valid sequence in the VCF file
    inline size_t getNumVariantsFromIndex(const string &vcffilename);

    void writeBCFRecords(vector<tbb::concurrent_bounded_queue<bcf1_t*>> &recq, htsFile *bcfout, bcf_hdr_t *bcfhdr, size_t numrecs);

    int command_argc;
    char** command_argv;

    MapInterpolater mapint;

    const string refFileName;
    vector<BooleanVector> reference; // reference haplotypes on target sites (false = ref, true = alt), size is NrefhapsxM (capacity: NrefhapsmaxxM), sample-major (includes space for phased targets for multiple phasing iterations)
    BooleanVector::data_type *refdata = 0;
    RingBuffer<BooleanVector> referenceT; // reference haplotypes on target sites (transposed), size MxNrefhaps (capacity: MxNrefhapsmax), haploids are encoded homozygous diploid!
    BooleanVector::data_type *refdataT = 0;


    // ATTENTION! Be aware that the underlying data chunks are not a combined portion of memory, but individually for each variant!
    // They also have to be freed individually
    vector<BooleanVector> referenceFullT; // complete reference haplotypes in chunk (false = ref, true = alt), size is Mrefx2xNref, SNP-major (no space for phased targets!)
    // capacity for haps in one sample
    size_t hapcapacity;

    vector<GenotypeVector> targets; // target genotypes (0 = hom ref, 1 = het, 2 = hom alt, 9 = missing), size is Ntarget
    vector<RingBuffer<size_t>> tgtmiss; // positions of missings in targets
    vector<RingBufferBool> tgtinphase; // input target phases, true == mat:1 pat:0, false others (only stored if --skipPhasing is enabled)
    vector<uint64_t> chrBpsReg; // base pair positions of target SNPs on target chromosome, size is Mreg
    RingBuffer<fp_type> cMs; // size is M
    fp_type cM0; // for the statistics: store first cM value
    vector<string> targetIDs; // identifier strings for target samples
    RingBufferBool isPhased; // size corresponds to the number of SNPs which will be written to phased output file in curr chunk (M if "!outputUnphased") and indicates if the entry is/will be phased or not
    vector<bool> haploidsTgt; // flag for all samples in target that are haploid
    vector<bool> haploidsTgt_initialized; // flag for all samples if the information in haploidsTgt is valid, usually set after first target SNP (just in the case if the first genotype is missing, it will be set later)
    vector<size_t> haploidsTgtMap; // map index of any haplotype vector (indexed from where haploids are encoded as haploids) to index in target haps where haploids are encoded as homozygous diploid
    size_t       nHaploidsTgt = 0; // number of haploid target samples
    vector<bool> haploidsRef; // flag for all samples in reference that are haploid (only allowed (and used) for X and Y chromosome)
    vector<bool> haploidsRef_initialized; // flag for all samples if the information in haploidsRef is valid, usually set after first reference SNP (just in the case if the first genotype is missing, it will be set later)
    vector<size_t> haploidsRefMap; // map index of any haplotype vector (indexed from where haploids are encoded as haploids) to index in reference where haploids are encoded as homozygous diploid
    size_t       nHaploidsRef = 0; // number of haploid reference samples
    RingBuffer<float> alleleFreqsCommon; // size corresponds to phased sites, contains the allele frequencies of each variant (-1.0 if unknown)
    RingBuffer<bcf1_t*> bcf_pout; // all records w/o genotypes for phased output; if "outputUnphased" is set, this vector keeps the records of the unphased SNPs including genotypes as well; ensure to add copies only!
    bcf_hdr_t *tgt_hdr_cpy; // copy of the target header used for the phased output file -> for imputation, the sample names are copied

    RingBufferBool isImputed; // size corresponds to imputation output of current chunk, indicates if the SNP is/will be imputed or not
    // for each analyzed target SNP the index in the full reference vector (relative to current chunk) that corresponds to this SNP, size is M
    RingBuffer<size_t> indexToRefFull;
    // for each reference SNP the index for the next common site in the the target (isPhased or bcf_pout), size as Mref
    // (index points to position in bcf_pout vector;
    //  if the current site is common, the index corresponds to the same site;
    //  if no common index follows, the stored value is M)
    RingBuffer<size_t> nextPidx;
    // for each reference SNP the index for the commonly analyzed target site (i.e. the index in phasedTarget), it points to the next common site if this site was not analyzed in phasing
    RingBuffer<size_t> ptgtIdx;

    vector<int64_t> positionsFullRefRegion; // chromosome positions of complete reference in region, 0-based
    vector<float> alleleFreqsFullRefRegion; // allele frequencies of each variant in complete reference in region
    vector<string> allelesFullRefRegion; // alleles, alternating ref and alt allele per variant, of complete reference in region
    vector<string> variantIDsFullRefRegion; // variant IDs of complete reference in region
    BooleanVector multiAllFlagsFullRefRegion;  // marked true if multi-allelic variant was split to bi-allelic, of complete reference in region

    size_t Nref;
    size_t Nrefhaps; // number of reference haplotypes (NOTE: originally haploids are encoded homozygous diploid!, so this is 2*Nref)
    size_t Ntarget; // number of target samples
    size_t Nrefhapsmax; // maximum number of reference samples (Nrefhaps if iters <= 1, Nrefhaps+Ntargethaps else)
    size_t M; // number of SNPs common in target and reference in current chunk (incl. overlap)
    size_t Mglob; // complete number of SNPs in target file
    size_t Mref; // number of SNPs in reference (for imputation) in current chunk (incl. overlap)
    size_t Mrefreg; // number of SNPs in reference (for imputation) in selected region -- only set with a Qref!
    size_t Mrefglob; // global number of SNPs in reference (before reduction to region)
    size_t MrefMultiAllreg; // number of SNPs in reference that originate from a split multi-allelic variant in selected region -- only set with a Qref!
    size_t MrefMultiAllglob; // global number of SNPs in reference that originate from a split multi-allelic variant (before reduction to region) -- only set with a Qref!
    size_t MLeftOv, MRightOv; // number of common tgt sites in overlap to previous (left) or next (right) chunk
    size_t MrefLeftOv, MrefRightOv; // number of ref sites in overlap to previous (left) or next (right) chunk
    size_t bunchsize; // number of sites for each file to be imputed in parallel (nbunches * bunchsize * num_files = Mref)
    size_t nbunches = 1;
    bool useExclude; // using an exclude file?

    // continuously updated while reading, after reading used to update M/Mref: currM: number of common tgt sites in chunk (incl. overlap); currMref: number of reference sites in chunk (incl. overlap)
    size_t currM = 0, currMref = 0;
    // after reading: number of tgt sites in chunk in overlap to NEXT chunk.
    size_t currMOverlap = 0;

    bool allowRefAltSwap;
    bool allowStrandFlip;
    bool excludeMultiAllRef;
    int chrom;
    string chromlit;
    bool chrliterally;
    bool setChrom = false;
    bool checkedChrom = false;
    // global region (end inclusive!), 1-based!
    int64_t startRegionBp, endRegionBp;
    int64_t startFlankRegionBp, endFlankRegionBp; // including flanks
    // this is the maximum size of a chunk in bytes
    const uint64_t maxchunkmem;
    const uint64_t chunkflanksize; // number of tgt variants in the overlap to the left and to the right (total overlap = 2xchunkflanksize)
    // region for current chunk (without flanks) (end inclusive!), 0-based!
    // set after reading the corresponding target data
    int64_t startChunkBp, endChunkBp;
    vector<size_t> endChunkFlankIdx; // tgt index for end of chunk (with flanking region, exclusive)
    // for phasing: start phase for each target
    // set for next chunk while writing phased or imputed data according to first heterozygous site in overlap region
    vector<bool> chunkStartPhase;

    bool outputUnphased; // if phased output is written, also write the unphased sites?
    string outFileName;
    string outFileNameConf;
    string outFileNameImp;
    string outFileNameInfo;
    string writeMode;

    uint32_t iters; // number of phasing iterations

    bool doImputation; // is imputation also done? -> keep complete reference etc
    bool noImpMissing; // will missings be imputed during phasing?
    bool skipPhasing;

    bool noMissingIDs; // convert missing variant IDs to chr:pos:ref:alt?

    // random number generator
    default_random_engine randGen;
    uniform_real_distribution<double> randDist;
    bool deterministic;

    // for imputation
    bool writeADosage = false;
    bool writeGDosage = false;
    bool writeProbs = false;
    float impR2filter = 0.0;
    float impMAFfilter = 0.0;
    unsigned num_workers = 1;
    unsigned num_files = 1;
    vector<size_t> num_sites_per_file;
    vector<size_t> site_offsets;
    vector<vector<float>> ads_vec;
    vector<vector<float>> ds_vec;
    vector<vector<float>> gp_vec;
    vector<int> mtgt_gt_vec;
    vector<void*> tgt_gt_vec;
    vector<vector<tbb::concurrent_bounded_queue<bcf1_t*>>> recqs;
    vector<thread> wrts;
    vector<htsFile*> bcfouts;
    bcf_hdr_t *imp_hdr; // header used for imputation output

    bool loadQuickRef;
    bool createQRef;
    string qreffilename;  // will be initialized only if we create a Qref
    ofstream qrefvarsofs; // will be initialized only if we create a Qref

    bool overwriteCalls;
    bool improveCalls;

    bool skipHeader;
//    int overrideChunks;

    size_t tgtlinesread;
    size_t reflinesread;
    size_t qrefregfirstidx = 0;
    size_t qrefreglastidx = UINT64_MAX; // exclusive!
    size_t qrefcurridxglob = 0;
    size_t qrefcurridxreg = 0; // points to the next variant to be loaded (equals the number of loaded variants in the region)
    ifstream qin; // for Qref
    bcf_srs_t *sr; // synchronized reader, for VCF/BCF
    bcf_hdr_t *ref_hdr; // header for reading reference VCF
    bcf_hdr_t *tgt_hdr; // header for reading target VCF
    int nChunks = 1; // total number of chunks (immediate number, only valid after all data was read)
    int minNChunks = 1; // minimum number of chunks the analysis is conducted with
    size_t maxChunkTgtVars = ~0ull; // maximum number of tgt variants allowed in a chunk
    int currChunk = -1; // indicates the currently loaded chunk, -1 indicates "chunk not yet loaded"
    size_t currChunkOffset = 0; // index offset of current chunk's first local imputation index to reference index in current region

    ofstream infofile; // imputation info file that contains information on each target variant

    size_t numThreads;
    bool usefpga;

    size_t conflictcnt = 0;

    size_t improvecntnochange = 0;
    size_t improvecnthet2het = 0;
    size_t improvecnthom2hom = 0;
    size_t improvecnthet2hom = 0;
    size_t improvecnthom2het = 0;
    double totalimprovenochange = 0.0;
    double totalimprovehet2het = 0.0;
    double totalimprovehom2hom = 0.0;
    double totalimprovehet2hom = 0.0;
    double totalimprovehom2het = 0.0;

    struct VCFStats {
        size_t M = 0, Mref = 0;
        size_t MrefMultiAllreg = 0;
        size_t Mexclude = 0, MtargetOnly = 0, MrefOnly = 0;
        size_t MmultiAllelicTgt = 0, MmultiAllelicRefTgt = 0, MmultiAllelicRefOnly = 0, MmultiSplittedTgt = 0, MmultiSplittedRefOnly = 0, MmonomorphicRefTgt = 0;
        size_t MwithMissingRef = 0, MwithUnphasedRef = 0, Munphased = 0;
        size_t MrefAltError = 0, numRefAltSwaps = 0, numStrandFlips = 0, numRefAltSwapAndStrandFlip = 0;
        size_t GmissingRef = 0, GunphasedRef = 0, GmissingTarget = 0, GmultiFilledTgt = 0, GmultiFilledRefOnly = 0;
        size_t MmissingRefOnly = 0, GmissingRefOnly = 0, MunphasedRefOnly = 0, GunphasedRefOnly = 0;
        VCFStats& operator+=(const VCFStats& s) {
            M += s.M; Mref += s.Mref;
            MrefMultiAllreg += s.MrefMultiAllreg;
            Mexclude += s.Mexclude; MtargetOnly += s.MtargetOnly; MrefOnly += s.MrefOnly;
            MmultiAllelicTgt += s.MmultiAllelicTgt; MmultiAllelicRefTgt += s.MmultiAllelicRefTgt; MmultiAllelicRefOnly += s.MmultiAllelicRefOnly;
            MmultiSplittedTgt += s.MmultiSplittedTgt; MmultiSplittedRefOnly += s.MmultiSplittedRefOnly; MmonomorphicRefTgt += s.MmonomorphicRefTgt;
            MwithMissingRef += s.MwithMissingRef; MwithUnphasedRef += s.MwithUnphasedRef; Munphased += s.Munphased;
            MrefAltError += s.MrefAltError; numRefAltSwaps += s.numRefAltSwaps; numStrandFlips += s.numStrandFlips; numRefAltSwapAndStrandFlip += s.numRefAltSwapAndStrandFlip;
            GmissingRef += s.GmissingRef; GunphasedRef += s.GunphasedRef; GmissingTarget += s.GmissingTarget; GmultiFilledTgt += s.GmultiFilledTgt; GmultiFilledRefOnly += s.GmultiFilledRefOnly;
            MmissingRefOnly += s.MmissingRefOnly; GmissingRefOnly += s.GmissingRefOnly; MunphasedRefOnly += s.MunphasedRefOnly; GunphasedRefOnly += s.GunphasedRefOnly;
            return *this;
        }
    };
    VCFStats globalstats;

    bool yaml = false;

    static const int QREFV_MAJ = 1;
    static const int QREFV_MIN = 1;

    static const int TMPBUFSIZE = 128;

    static const int CHRX = 23; // mapping chromosome X to number 23
    static const int CHRY = 24; // mapping chromosome Y to number 24

};

#endif /* VCFDATA_H_ */
