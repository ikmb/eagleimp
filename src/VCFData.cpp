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

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <tbb/concurrent_queue.h>
#include <thread>
#include <functional>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <random>

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <version.h>

#include "MapInterpolater.h"
#include "StatusFile.h"
#include "utils.h"
#include "Stopwatch.h"

#include "VCFData.h"
#include "hybridsys/ThreadUtils.h"

using namespace placeholders;

VCFData::VCFData(const Args& args, int argc, char** argv, const vector<FPGAConfigurationEagleImp> &fpgaconfs)
    : command_argc(argc), command_argv(argv),
      refFileName(args.ref),
      allowRefAltSwap(args.allowRefAltSwap),
      allowStrandFlip(args.allowStrandFlip),
      excludeMultiAllRef(args.excludeMultiAllRef),
      maxchunkmem(args.maxChunkMemory * 1024*1024*1024), // value is provided in GiB from Args
      chunkflanksize(args.chunkFlankSize),
      outputUnphased(args.outputUnphased), // only true if phased output will be written as well
      writeMode(args.writeMode),
      iters(args.iters),
      doImputation(!args.skipImputation),
      noImpMissing(!args.impMissing),
      skipPhasing(args.skipPhasing),
      noMissingIDs(args.noMissingIDs),
      randGen(chrono::system_clock::now().time_since_epoch().count()), // random generator seed
      randDist(0.0, 1.0), // random generator distribution
      deterministic(args.deterministic),
      writeADosage(args.writeADosage),
      writeGDosage(args.writeGDosage),
      writeProbs(args.writeProb),
      impR2filter(args.impR2filter),
      impMAFfilter(args.impMAFfilter),
      loadQuickRef(args.isQRef),
      createQRef(args.createQuickRef),
      overwriteCalls(args.overwriteCalls),
      improveCalls(args.improveCalls),
      skipHeader(args.skipHeader),
//      overrideChunks(args.overrideChunks),
      numThreads(args.num_threads),
      usefpga(!fpgaconfs.empty()),
      yaml(args.yaml)
      {

    startRegionBp = args.bpStart;
    endRegionBp = args.bpEnd;
    startFlankRegionBp = args.bpStart <= 1+args.bpFlanking ? 1 : args.bpStart - args.bpFlanking;
    endFlankRegionBp = MAX_CSI_COOR-1 - args.bpEnd <= args.bpFlanking ? MAX_CSI_COOR-1 : args.bpEnd + args.bpFlanking;

    // process metadata of target and reference files
    // target will be ignored when Qref is created
    processMeta(refFileName, args.vcfTarget, args.vcfExclude, args, fpgaconfs);

    if (!createQRef) {

        outFileName = args.outPrefix;
        outFileName.append(".phased");
        outFileName.append(getOutputSuffix());
        outFileNameConf = args.outPrefix;
        outFileNameConf.append(".phased.confidences");
        outFileNameImp = args.outPrefix;
        outFileNameImp.append(".imputed");
        outFileNameImp.append(getOutputSuffix());
        outFileNameInfo = args.outPrefix;
        outFileNameInfo.append(".varinfo");

        // read genetic map (and test for errors)
        if (!skipPhasing) {
//            if (chrom != CHRY && args.genMapFile.empty()) {
            if (args.genMapFile.empty()) { // if you provide a haploid file, either specify a genetic map (although not needed) or explicit state --skipPhasing
                StatusFile::addError("No genetic map file specified. Try --help.");
                exit(EXIT_FAILURE);
            }
            mapint = MapInterpolater(args.genMapFile, chrom); // chromosomes with literal-only identifiers need a genetic map without "chr" column or 0 as integer identifier
        }
    } else { // createQRef
        // qref filename
        if (refFileName.size() >= 7 && refFileName.substr(refFileName.size()-7).compare(".vcf.gz") == 0) // special case for .vcf.gz
            qreffilename = refFileName.substr(0, refFileName.size() - 7);
        else // standard case: simply leave out file suffix starting with '.'
            qreffilename = refFileName.substr(0, refFileName.rfind('.'));
        qreffilename.append(".qref");

        // open output stream for variants
        qrefvarsofs.open(qreffilename+".vars", ios_base::binary);
    }

}

inline void VCFData::processMeta(const string &refFile, const string &vcfTarget, const string &vcfExclude, const Args &args, const vector<FPGAConfigurationEagleImp> &fpgaconfs) {

    M = 0;
    Mglob = 0;
    Mref = 0;
    Mrefreg = 0;
    Mrefglob = 0;
    MrefMultiAllreg = 0;
    MrefMultiAllglob = 0;
    MLeftOv = 0;
    MRightOv = 0;
    MrefLeftOv = 0;
    MrefRightOv = 0;

    cout << "\nAnalysing input files..." << endl;

    useExclude = !createQRef && !vcfExclude.empty();
    if (useExclude)
        cout << "  Using VCF exclusion file." << endl;

    cout << endl;

    // prepare status
    StatusFile::addInfo("<h3>General information:</h3>", false);
    // only for YAML info
    stringstream yamlinfo;
    yamlinfo << "\n";

    // determine number of SNPs in input file from corresponding index
    if (!createQRef) {
        Mglob = getNumVariantsFromIndex(vcfTarget);
        chrBpsReg.reserve(Mglob);
    }

    // We learn the chromosome from the reference file name per convention, i.e. it has to start with the number or with "chr" followed by the number,
    // and the identifier stops with a '.' or '_'. Completely literal names are also allowed now and will be referred to as number 0.
    // TODO probably we do not need this convention in the future?
    {
		string chrname(refFile);
		size_t p = refFile.rfind('/');
		if (p != string::npos) // need to cut off the path
			chrname = refFile.substr(p+1);
		p = chrname.find_first_of("._"); // find the stop character
		chrname = chrname.substr(0, p);
		chrom = chrNameToNumber(chrname, chromlit, &chrliterally);
		// if we load a Qref, the number and literal name is fixed here.
		// with a VCF reference, the name and number are checked to be consistent with the filename,
		// but the actual contents are taken later from inside the file
		setChrom = loadQuickRef;
    }

    // initialize VCF reader
    sr = bcf_sr_init();

//    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, allowRefAltSwap ? BCF_SR_PAIR_BOTH_REF : BCF_SR_PAIR_EXACT);
    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF);
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    bcf_sr_set_threads(sr, numThreads);

    {
        // set the region (including flanks)
        stringstream s;
        if (chrom > 0) { // we have a valid chromosome number (otherwise, we only have a literal identifier)
			s << chrom << ":" << startFlankRegionBp << "-" << endFlankRegionBp;
			s << ",chr" << chrom << ":" << startFlankRegionBp << "-" << endFlankRegionBp;
			if (chrom == CHRX) { // we literally add X as well
				s << ",X:" << startFlankRegionBp << "-" << endFlankRegionBp;
				s << ",chrX:" << startFlankRegionBp << "-" << endFlankRegionBp;
			}
			if (chrom == CHRY) { // or Y...
				s << ",Y:" << startFlankRegionBp << "-" << endFlankRegionBp;
				s << ",chrY:" << startFlankRegionBp << "-" << endFlankRegionBp;
			}
        } else { // only literal identifier
        	s << chromlit << ":" << startFlankRegionBp << "-" << endFlankRegionBp;
        	s << ",chr" << chromlit << ":" << startFlankRegionBp << "-" << endFlankRegionBp;
        }
        if (bcf_sr_set_regions(sr, s.str().c_str(), 0) != 0) {
            StatusFile::addError("Failed to initialize the region: " + s.str());
            exit(EXIT_FAILURE);
        }

    }

    if (!loadQuickRef && !bcf_sr_add_reader(sr, refFile.c_str())) {
        StatusFile::addError("Could not open reference file for reading: " + string(bcf_sr_strerror(sr->errnum)));
        exit(EXIT_FAILURE);
    }

    if (!createQRef) {
        if (!bcf_sr_add_reader(sr, vcfTarget.c_str())) {
            StatusFile::addError("Could not open target file for reading: " + string(bcf_sr_strerror(sr->errnum)));
            exit(EXIT_FAILURE);
        }
    }

    if (useExclude && !bcf_sr_add_reader(sr, vcfExclude.c_str())) {
      StatusFile::addError("Could not open exclude file for reading: " + string(bcf_sr_strerror(sr->errnum)));
      exit(EXIT_FAILURE);
    }

    tgt_hdr = NULL;
    tgt_hdr_cpy = NULL;
    if (!createQRef) {
        tgt_hdr = bcf_sr_get_header(sr, loadQuickRef ? 0 : 1);
        tgt_hdr_cpy = bcf_hdr_dup(tgt_hdr); // copy and keep for writing later
        appendVersionToBCFHeader(tgt_hdr_cpy);
        Ntarget = bcf_hdr_nsamples(tgt_hdr);

        haploidsTgt.clear();
        haploidsTgt.resize(Ntarget, false); // initialized with all diploid
        haploidsTgt_initialized.clear();
        haploidsTgt_initialized.resize(Ntarget, false); // not initialized yet
        haploidsTgtMap.clear();
        haploidsTgtMap.reserve(2*Ntarget);

        // Read target sample IDs
        targetIDs.reserve(Ntarget);
        for (size_t i = 0; i < Ntarget; i++)
            targetIDs.push_back(tgt_hdr->samples[i]);

        // framework for target data
        targets.clear();
        targets.resize(Ntarget);
        tgtmiss.clear();
        tgtmiss.resize(Ntarget);
        tgtinphase.clear();
        tgtinphase.resize(Ntarget);

        chunkStartPhase.clear();
        chunkStartPhase.resize(Ntarget, false); // w.l.o.g. the first starting phase is the reference allele for the maternal path
    }

    ref_hdr = NULL;
    if (loadQuickRef) { // load metadata from Qref
        qRefOpenReadMeta();
    } else { // prepare reading from VCF reference for creating a Qref
        ref_hdr = bcf_sr_get_header(sr, 0);
        Nref = bcf_hdr_nsamples(ref_hdr);
        Nrefhaps = 2*Nref;
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
        Nrefhapsmax = iters > 1 ? Nrefhaps + 2*(DEBUG_TARGET_STOP-DEBUG_TARGET_START+1) : Nrefhaps;
#else
        Nrefhapsmax = iters > 1 ? Nrefhaps + 2*Ntarget : Nrefhaps;
#endif
        hapcapacity = roundToMultiple<size_t>(Nrefhaps, UNITWORDS*sizeof(BooleanVector::data_type)*8)/8; // space for 2*Nref haps

        haploidsRef.clear();
        haploidsRef.resize(Nrefhapsmax/2, false); // initialized with all diploid
        haploidsRef_initialized.clear();
        haploidsRef_initialized.resize(Nref, false); // not initialized yet
        haploidsRefMap.clear();
        haploidsRefMap.reserve(Nrefhapsmax);

        // determine number of SNPs in input file by reading the index
        Mrefglob = getNumVariantsFromIndex(refFile);

        // reserve space for globals
        size_t Mrefpre = Mrefglob + Mrefglob/10; // 10% extra for multi-allelic splits
        alleleFreqsFullRefRegion.reserve(Mrefpre); // no effect if already loaded
        positionsFullRefRegion.reserve(Mrefpre);
        allelesFullRefRegion.reserve(2*Mrefpre);
        variantIDsFullRefRegion.reserve(Mrefpre);
        size_t macapacity = roundToMultiple<size_t>(Mrefpre, UNITWORDS * sizeof(BooleanVector::data_type)*8)/8; // space for Mrefpre SNPs
        BooleanVector::data_type* madata = (BooleanVector::data_type*) MyMalloc::malloc(macapacity, string("madata_processMeta"));
        multiAllFlagsFullRefRegion.setDataAndInit(madata, macapacity, 0, false);
    }

    stringstream s;
    s << "<table><tr><td>Chromosome:         </td><td>" << chromlit << "</td></tr>";
    StatusFile::addInfo(s.str(), false);
    yamlinfo << "    Chromosome: " << chromlit << "\n";

    if (!createQRef) {
        StatusFile::addInfo("<tr><td>Target samples:     </td><td>" + to_string(Ntarget) + "</td></tr>", false);
        StatusFile::addInfo("<tr><td>Target variants:    </td><td>" + to_string(Mglob) + "</td></tr>", false);
        yamlinfo << "    Target samples: " << Ntarget << "\n";
        yamlinfo << "    Target variants: " << Mglob << "\n";
    }
    StatusFile::addInfo("<tr><td>Reference samples:  </td><td>" + to_string(Nref) + "</td></tr>", false);
    yamlinfo << "    Reference samples: " << Nref << "\n";

    if (Mrefreg == Mrefglob || !loadQuickRef) {
        StatusFile::addInfo("<tr><td>Reference variants: </td><td>" + to_string(Mrefglob) + "</td></tr>", false);
        yamlinfo << "    Reference variants: " << Mrefglob << "\n";
    } else {
        StatusFile::addInfo("<tr><td>Reference variants (global): </td><td>" + to_string(Mrefglob) + "</td></tr>", false);
        StatusFile::addInfo("<tr><td>Reference variants (region): </td><td>" + to_string(Mrefreg) + "</td></tr>", false);
        yamlinfo << "    Reference variants (global): " << Mrefglob << "\n";
        yamlinfo << "    Reference variants (region): " << Mrefreg << "\n";
    }
    StatusFile::addInfo("</table>", false);
    if (yaml)
        StatusFile::addInfoYAML("General information", yamlinfo.str());

    if (Mrefreg <= 1ull && loadQuickRef) { // no data in reference
        StatusFile::addWarning("<b>Analysis skipped:</b> Too few reference data for this region. (Mrefreg = " + to_string(Mrefreg) + ")");
        StatusFile::updateStatus(1, "Skipped");
        exit(EXIT_SUCCESS);
    }

    // Memory requirements check:
    // We estimate the amount of required memory and perform a division into chunks if necessary and possible.
    // We also check, if the data can be used with the current FPGA configuration (if available)
    {
        // FPGA check:
        if (!createQRef && usefpga) { // we have an FPGA configuration
            // take only the first conf
            const auto &conf = fpgaconfs[0];

            // check maximum haps
            if (conf.getMaxHaps() < Nrefhapsmax) {
                stringstream ss;
                if (args.iters > 1 && conf.getMaxHaps() >= Nrefhaps) {
                    ss << "<b>Disabling FPGA:</b> Analysis requires support for " << Nrefhapsmax << " haplotypes due to several phasing iterations,<br>\n"
                       << "  but FPGA only supports " << conf.getMaxHaps() << ". To use the FPGA try only one phasing iteration by applying -i 1.";
                } else {
                    ss << "<b>Disabling FPGA:</b> Analysis requires support for " << Nrefhapsmax << " haplotypes,<br>\n"
                       << "  but FPGA only supports " << conf.getMaxHaps() << ".";
                }
                StatusFile::addWarning(ss.str());
                usefpga = false;
            }

            // check maximum K
            if (conf.getMaxK() < min(args.K, Nrefhapsmax)) {
                stringstream ss;
                if (args.iters > 1 && conf.getMaxK() >= min(args.K, Nrefhaps)) {
                    ss << "<b>Disabling FPGA:</b> Analysis requires support for K=" << min(args.K, Nrefhapsmax) << " conditioning haplotypes<br>\n"
                       << "  due to several phasing iterations, but FPGA only supports " << conf.getMaxK() << ".<br>\n"
                       << "  Try only one phasing iteration by applying -i 1 or a lower K by -K " << conf.getMaxK() << ".";
                } else {
                    ss << "<b>Disabling FPGA:</b> Analysis requires support for K=" << min(args.K, Nrefhapsmax) << " conditioning haplotypes,<br>\n"
                       << "  but FPGA only supports " << conf.getMaxK() << ". Try applying a lower K, e.g. -K " << conf.getMaxK() << ".";
                }
                StatusFile::addWarning(ss.str());
                usefpga = false;
            }
        }

        // condensed references + PBWTs (estimated, if every third site is a call site, fwd+bck and ref+inc, no cref required when using FPGA,
        // but number corresponds to all FPGA pipelines times the capacity of the FPGA processor outqueue, which is fixed to 2)
        size_t reqphasecref = (skipPhasing ? 0 : (usefpga ? fpgaconfs[0].getNumPipelines() * 2 : numThreads)) * min(Nrefhapsmax, args.K) * (Mglob/3) / (usefpga ? 2 : 1);
        size_t reqphasedos = Ntarget * Mglob * 8; // phased dosages
        size_t reqimpref = doImputation ? (Nrefhaps * Mrefglob) / 8 : 0ull; // required size in bytes for reference haplotypes
        size_t reqimppbwt = doImputation ? Nrefhaps * Mglob : 0ull; // PBWT and refT of common reference (if every tgt site is also found in the reference)
        size_t reqimpabsp = doImputation ? Nrefhaps * Mglob * 4 : 0ull; // absolute permutation arrays for reference PBWT (if every tgt site is also found in the reference)

        size_t reqimpbunchhaps = 0ull;
        size_t reqimpqueue = 0ull;

        if (doImputation) {
            num_files = args.num_files;
            num_workers = max(1u, args.num_threads - num_files);

            // prepare number of sites per temp file for imputation output
            // - first and last field represent the vars in the overlap that should be ignored
            // - will be filled after chunk processing has finished, but first field stays zero for the first chunk
            num_sites_per_file.clear();
            num_sites_per_file.resize(num_files+2);

            // (preliminary) bunchsize in number of variants:
            // estimated from user parameter and memory requirements, assumption based on equally distributed variants;
            // will be adapted to real chunk size later
            bunchsize = max((size_t)1, (args.ibunchsize*(size_t)1000) / Ntarget); // minimum bunch size is 1
            bunchsize = roundToMultiple(bunchsize, (size_t)num_workers); // important to have an equal load on all worker threads

            reqimpbunchhaps = Ntarget * bunchsize / 4; // imputation haplotypes (mat+pat) in a bunch (in bytes)
            reqimpqueue = Ntarget * bunchsize * num_files * 2; // number of records in output queues x number of target samples per record
            size_t entrysize = 8; // size of two ints for the haplotypes
            if (writeADosage)
                entrysize += 8; // size of two floats for the allele dosages
            if (writeGDosage)
                entrysize += 4; // size of one float for the genotype dosage
            if (writeProbs)
                entrysize += 12; // size of three floats for each genotype probability
            reqimpqueue *= entrysize;

            // check if required static size does not exceed half of the chunk memory, otherwise reduce bunch size
            if (reqimpbunchhaps + reqimpqueue > maxchunkmem/2) {
                bunchsize /= divideRounded(reqimpbunchhaps + reqimpqueue, maxchunkmem/2);
                if (bunchsize < num_workers) { // we need at least one record per bunch per worker!
                    StatusFile::addError("Too many samples. Sorry.");
                    exit(EXIT_FAILURE);
                }
                // only set to multiple of num_workers if the memory increase would not be too much!
                if (bunchsize > 4*num_workers)
                    bunchsize = roundToMultiple(bunchsize, (size_t)num_workers);
                // recalculate required sizes
                reqimpbunchhaps = Ntarget * bunchsize / 4; // imputation haplotypes (mat+pat) in a bunch (in bytes)
                reqimpqueue = Ntarget * bunchsize * num_files * 2 * entrysize; // number of records in output queues x number of target samples per record x size of the entry
            }
        }
        size_t reqsum_dyn = reqphasecref + reqphasedos + reqimpref + reqimppbwt + reqimpabsp;
        size_t reqsum_stat = reqimpbunchhaps + reqimpqueue;
        size_t reqsafety = reqsum_dyn/10;
        reqsum_dyn += reqsafety;

        if (!createQRef) { // give a short memory summary
//            if (args.debug) {
                cout << "Roughly estimated maximum memory requirements for this analysis (without region settings, if every 3rd site is a call site):" << endl;
                cout << "  Chunk memory:" << endl;
                cout << "    Crefs + PBWTs  : " << divideRounded(reqphasecref, 1024ul*1024ul) << " MiB" << endl;
                cout << "    Phased dosages : " << divideRounded(reqphasedos, 1024ul*1024ul) << " MiB" << endl;
                cout << "    Reference haps : " << divideRounded(reqimpref, 1024ul*1024ul) << " MiB" << endl;
                cout << "    RefT + PBWT    : " << divideRounded(reqimppbwt, 1024ul*1024ul) << " MiB" << endl;
                cout << "    PBWT absPerm   : " << divideRounded(reqimpabsp, 1024ul*1024ul) << " MiB" << endl;
                cout << "    10% safety     : " << divideRounded(reqsafety, 1024ul*1024ul) << " MiB" << endl;
                cout << "  Static (chunk independent) memory:" << endl;
                cout << "    Imputed haps   : " << divideRounded(reqimpbunchhaps, 1024ul*1024ul) << " MiB" << endl;
                cout << "    Outqueue space : " << divideRounded(reqimpqueue, 1024ul*1024ul) << " MiB" << endl;
                cout << "  SUM              : " << divideRounded(reqsum_dyn + reqsum_stat, 1024ul*1024ul) << " MiB" << endl;
//            } else
//                cout << "Roughly estimated maximum memory requirements for this analysis: " << divideRounded(reqsum_dyn + reqsum_stat, 1024ul*1024ul) << " MiB" << endl;
            cout << "Maximum allowed memory: " << divideRounded(maxchunkmem, 1024ul*1024ul) << " MiB" << endl;
        }


        // Determine number of chunks: TODO check influence of region settings
        // If the user's set a region, we don't know how many variants are actually in that region,
        // so we generate chunk sizes as if no region was set, but reading the region will
        // start filling the first chunk and so on with risking the last chunks to be empty.
        // The last chunk with data will most likely not be filled completely, but the overlap
        // will ensure that data processing can still be done.

        if (createQRef) { // no chunking required for creating a Qref
            minNChunks = 1;
            nChunks = 1;
        } else {
//            if (overrideChunks) { // override chunk calculation
//                nChunks = overrideChunks;
//                StatusFile::addWarning("Automatic chunk division disabled by user.");
//                if (reqsum_dyn/nChunks > maxchunkmem-reqsum_stat) {
//                    StatusFile::addWarning("The analysis of your target data together with the reference will probably not fit into memory!<br>\n"
//                            "  Estimated maximum memory requirements per chunk: " + to_string(reqsum_dyn/nChunks+reqsum_stat) + " bytes.");
//                }
//            } else {

//                nChunks = divideRounded(reqsum_dyn, maxchunkmem-reqsum_stat); // a safety margin is already included in reqsum
//                if (nChunks > 1) {
//                    // recalculate memory requirements per chunk including overlap
//                    size_t reqsum_dyn_chunk = (Mglob/nChunks + 2*chunkflanksize) * reqsum_dyn / Mglob; // estimated from the relation Mchunk/Mglob
//                    while (reqsum_dyn_chunk > maxchunkmem-reqsum_stat && Mglob/nChunks >= chunkflanksize) {
//                        int nChunks_tmp = nChunks;
//                        nChunks = divideRounded(reqsum_dyn_chunk*nChunks, maxchunkmem-reqsum_stat);
//                        if (nChunks == nChunks_tmp) // to ensure an exit of the while loop
//                            nChunks++;
//                        reqsum_dyn_chunk = (Mglob/nChunks + 2*chunkflanksize) * reqsum_dyn / Mglob;
//                    }
//
//                }
            // The number of tgt sites that are allowed in a chunk at max:
            // required dynamic memory per site (assuming an equal distribution)
            size_t reqpersite = reqsum_dyn/Mglob;
            // additionally required dynamic memory per site if an overlap is required (more than one chunk, space for ref and refT on common sites is reserverd twice)
            size_t reqperovsite = reqimppbwt/Mglob;
            // max vars in a chunk adjusted for using overlaps and taking care of the required static memory
            maxChunkTgtVars = divideRounded(maxchunkmem-reqsum_stat, reqpersite+reqperovsite);

            cerr << "maxChunkTgtVars: " << maxChunkTgtVars << endl;

                // check if this estimation is conform with the FPGA configuration
                if (usefpga) {
                    // TODO need correction for FPGA?
//                    // determine number of chunks according to maximum supported sites
//                    int nChunks_old = nChunks;
//                    while (fpgaconfs[0].getMaxSites() < (Mglob/nChunks+2*chunkflanksize) && Mglob/nChunks >= chunkflanksize) { // the maximum number of supported sites is smaller than the sites in each chunk
//                        int nChunks_tmp = nChunks;
//                        nChunks = divideRounded((Mglob+2*nChunks*chunkflanksize), fpgaconfs[0].getMaxSites());
//                        if (nChunks == nChunks_tmp) // to ensure an exit of the while loop
//                            nChunks++;
//                    }
//                    if (nChunks != nChunks_old)
//                        StatusFile::addWarning("Increased number of chunks to " + to_string(nChunks) + " due to the FPGA's maximum number of sites restriction.");
//
//                    // check if FPGA memory is sufficient with the current number of chunks
//                    size_t maxKbytes = roundToMultiple(min(args.K, Nrefhapsmax), (size_t) 512) / 8;
//                    size_t maxMemPerSite = maxKbytes*4 * fpgaconfs[0].getNumPipelines(); // required memory on the FPGA if ~every 3rd site is a split site (fwd+bck and ref+inc per split site)
//                    size_t maxSites = (Mglob/nChunks + 2*chunkflanksize)/3; // if every third site is a split site
//                    int nChunks_old2 = nChunks;
//                    while (fpgaconfs[0].getAvailableRAM() < maxMemPerSite*maxSites && Mglob/nChunks >= chunkflanksize) {
//                        int nChunks_tmp = nChunks;
//                        nChunks = divideRounded(maxMemPerSite*(Mglob+2*nChunks*chunkflanksize)/3, fpgaconfs[0].getAvailableRAM());
//                        if (nChunks == nChunks_tmp) // to ensure an exit of the while loop
//                            nChunks++;
//                        maxSites = (Mglob/nChunks + 2*chunkflanksize)/3; // if every third site is a split site
//                    }
//                    if (nChunks != nChunks_old2)
//                        StatusFile::addWarning("Increased number of chunks to " + to_string(nChunks) + " due to the FPGA's available memory restriction.");
//
//                    // check, if we increased the number of chunks too much
//                    if (Mglob/nChunks < 2*chunkflanksize) {
//                        nChunks = nChunks_old;
//                        usefpga = false;
//                        StatusFile::addWarning("<b>Disabled FPGA:</b> Chunk size too small. Reverted number of chunks to " + to_string(nChunks) + " and continuing with CPU only.");
//                    }
                }
//            }

            // this is a rough estimation intended to be smaller than the final number of chunks!
            minNChunks = divideRounded(Mglob, maxChunkTgtVars); // no overlaps, equal distribution

            // check if chunk size is large enough:
            // as an overlap will be defined on both ends of a chunk, chunk size must be at least two overlaps!
            if (minNChunks > 1 && maxChunkTgtVars < 2*chunkflanksize) {
                string serr("<b>Analysis not possible:</b> Too many chunks.");
//                if (!overrideChunks)
                serr += " If possible, try larger chunk memory size with --maxChunkMemory.";
                StatusFile::addError(serr);
                exit(EXIT_FAILURE);
            }

            if (minNChunks > 1) {
                // we calculated without overlap above, but with more than one chunk we need to consider overlaps now
                int64_t remvars = Mglob - (minNChunks*maxChunkTgtVars - (minNChunks-1)*2*chunkflanksize);
                if(remvars > 0) {
                    minNChunks += divideRounded((size_t)remvars, maxChunkTgtVars-2*chunkflanksize);
                }
            }

            nChunks = minNChunks; // might be increased!
            // preliminary number of bunches
            nbunches = divideRounded(divideRounded(Mrefglob, (size_t)nChunks), bunchsize * num_files);
        }

        // chunks are generated based on the number of target variants
        startChunkIdx.clear();
        startChunkFlankIdx.clear();
        endChunkFlankIdx.clear();

        // we start with the first target index
//        size_t curridx = 0;
        startChunkIdx.push_back(0);
        startChunkFlankIdx.push_back(0);
        // we only prepare preliminary limits for the next chunk, others will be generated on-the-fly
        if (minNChunks == 1) { // memory estimation thinks we can handle the imputation without chunking -> and this calculation should be very resolute
            startChunkIdx.push_back(Mglob);
            startChunkFlankIdx.push_back(Mglob);
            endChunkFlankIdx.push_back(Mglob);
        } else {
            startChunkIdx.push_back(maxChunkTgtVars-chunkflanksize);
            startChunkFlankIdx.push_back(maxChunkTgtVars-2*chunkflanksize);
            endChunkFlankIdx.push_back(maxChunkTgtVars);
        }
//        for (int chunk = 1; chunk < nChunks; chunk++)
//        {
//            // equally distribute variants
//            curridx += Mglob/nChunks + ((size_t)chunk < Mglob%nChunks ? 1 : 0); // add size of current chunk to get start index of next chunk
//            startChunkIdx.push_back(curridx);
//            startChunkFlankIdx.push_back(curridx - chunkflanksize); // as the chunk size will be larger than the flank size, this index will exist
//            endChunkFlankIdx.push_back(curridx + chunkflanksize); // as the chunk size will be larger than the flank size, this index will exist as well
//        }
//        // calculate the maximum possible extension of a chunk
//        maxChunkExtension = 0; //Mglob - startChunkFlankIdx.back() - 2*chunkflanksize;
//        // sentinel at the end
//        startChunkIdx.push_back(Mglob);
//        startChunkFlankIdx.push_back(Mglob); // important to have no overlap here!
//        endChunkFlankIdx.push_back(Mglob);

        if (!createQRef) {
            StatusFile::addInfo(string("<p class='pinfo'>Data will be processed in ").append(nChunks > 1 ? "at least " : "").append(to_string(nChunks)).append(" chunk").append(nChunks > 1 ? "s." : " (no splitting required).</p>"), false);
            if (yaml) {
                StatusFile::addInfoYAML("Chunks", to_string(nChunks));
            }
//            // DEBUG
//            cout << "Chunks: " << endl;
//            for (int c = 0; c < nChunks; c++) {
//                cout << c << ": " << startChunkFlankIdx[c] << " " << startChunkIdx[c] << " " << endChunkFlankIdx[c] << endl;
//            }
//            // __DEBUG
        }

        // the chunk regions are set after loading each corresponding chunk
    }

    currChunk = -1;
    imp_hdr = NULL;
    tgtlinesread = 0;
    reflinesread = 0;
}


inline void VCFData::initImpHeader() {
    imp_hdr = bcf_hdr_init("w"); // creates header with the first two lines: VCF version + FILTER=...
    {
        // add the (required) chromosome info to the header
        // since we only impute chromosome-wise, the rid field in each bcf record (later) will always be 0 independent of the chromosome number
        stringstream contig;
        contig << "##contig=<ID=" << chromlit;
        // only known in advance if we loaded a Qref; TODO try to add this info (and the following) to the header at the end just before writing it
        contig << ",length=" << (loadQuickRef ? positionsFullRefRegion.back()+1 : MAX_CSI_COOR-1) << ">";
        bcf_hdr_append(imp_hdr, contig.str().c_str());
    }

    appendVersionToBCFHeader(imp_hdr);

    bcf_hdr_append(imp_hdr, "##INFO=<ID=RefPanelAF,Number=A,Type=Float,Description=\"Allele frequency in imputation reference panel\">");
    bcf_hdr_append(imp_hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency from allele dosages\">");
    bcf_hdr_append(imp_hdr, "##INFO=<ID=MAF,Number=A,Type=Float,Description=\"Estimated minor allele frequency from allele dosages\">");
    bcf_hdr_append(imp_hdr, "##INFO=<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy (R-square)\">");
    bcf_hdr_append(imp_hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
    bcf_hdr_append(imp_hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    bcf_hdr_append(imp_hdr, "##INFO=<ID=TYPED,Number=0,Type=Flag,Description=\"Site was genotyped and phased prior to imputation\">");
//    bcf_hdr_append(imp_hdr, "##INFO=<ID=TYPED_ONLY,Number=0,Type=Flag,Description=\"Site was genotyped but is not phased or imputed\">");
    bcf_hdr_append(imp_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

    if (writeADosage)
        bcf_hdr_append(imp_hdr, "##FORMAT=<ID=ADS,Number=.,Type=Float,Description=\"Allele dosage per haplotype (for phased or imputed sites only)\">");
    if (writeGDosage)
        bcf_hdr_append(imp_hdr, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage\">");
    if (writeProbs)
        bcf_hdr_append(imp_hdr, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">");

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
    string samplename("DBGTGT");
    for (int samp = DEBUG_TARGET_START; samp <= DEBUG_TARGET_STOP; samp++) {
        bcf_hdr_add_sample(imp_hdr, samplename.append(to_string(samp)).c_str());
    }
#else
    // take the same samples as in the original header
    int nsamp = bcf_hdr_nsamples(tgt_hdr);
    for (int s = 0; s < nsamp; s++)
        bcf_hdr_add_sample(imp_hdr, (tgt_hdr->samples)[s]);
#endif
}

// initialize variant info file
inline void VCFData::initInfoFile() {
    infofile.open(outFileNameInfo);
    // header
    infofile << "#CHR\tPOS\tID\tREF\tALT\tINCL/EXCL\tRID\tRREF\tRALT\tExplanation" << endl;
}

// add excluded variant to info file
inline void VCFData::addToInfoFileExcluded(bcf1_t* tgt, const string& explanation) {
    bcf_unpack(tgt, BCF_UN_STR); // for the ID and alleles
    infofile << chromlit << "\t" << tgt->pos+1 << "\t" << tgt->d.id << "\t" << tgt->d.allele[0] << "\t";
    if (tgt->n_allele > 1)
        infofile << tgt->d.allele[1];
    for (uint32_t i=2; i < tgt->n_allele; i++)
        infofile << "," << tgt->d.allele[i];
    infofile << "\tEXCLUDED\t\t\t\t" << explanation << endl;
}
inline void VCFData::addToInfoFileExcluded(size_t pos, const string& tgtID, const string& tref, const string& talt, const string& explanation) {
    infofile << chromlit << "\t" << pos+1 << "\t" << tgtID << "\t" << tref << "\t" << talt << "\tEXCLUDED\t\t\t\t" << explanation << endl;
}

// add included variant to info file
inline void VCFData::addToInfoFileIncluded(size_t pos, const string& tgtID, const string& tref, const string& talt, const string& refID, const string& rref, const string& ralt, const string& explanation) {
    infofile << chromlit << "\t" << pos+1 << "\t" << tgtID << "\t" << tref << "\t" << talt << "\tINCLUDED\t";
    if (!noMissingIDs || refID.compare("."))
        infofile << refID;
    else // the problem is, that the correction for the imputation will be done later, so we do it here again...
        infofile << chromlit << ":" << pos+1 << ":" << rref << ":" << ralt;
    infofile << "\t" << rref << "\t" << ralt << "\t" << explanation << endl;
}

// every data array re-allocated here is always addessed relative to the current chunk.
// exceptions are those arrays with metadata loaded from a Qref (which are marked as global)
void VCFData::processNextChunk() {

    currChunk++;

    // the statistics collected here are only for the current chunk and without the overlap from the previous chunk
    VCFStats lstats;

    if (!createQRef) {
        cout << "\n----------------------------------------------------------------\n--- ";
        stringstream s;
        s << "<h4>Chunk " << currChunk+1 << "</h4>";
        StatusFile::setContext(s.str());
        cout << "Chunk " << currChunk+1 << "/" << nChunks << ":" << endl;
        cout << "----------------------------------------------------------------" << endl;

        // DEBUG
        cerr << "sf/s/e: " << startChunkFlankIdx[currChunk] << "/" << startChunkIdx[currChunk] << "/" << endChunkFlankIdx[currChunk]
             << " sizef/size: " <<  endChunkFlankIdx[currChunk] - startChunkFlankIdx[currChunk] << "/" << endChunkFlankIdx[currChunk] - startChunkIdx[currChunk]
             << " tgtlinesread: " << tgtlinesread
             << " sfn/sn: " << startChunkFlankIdx[currChunk+1] << "/" << startChunkIdx[currChunk+1]
             << endl;

        //
        // prepare memory for target data of current chunk, copy overlap from previous chunk
        //

        // determine Mrefpre for chunks:
        // number of reference vars per chunk if tgt and ref vars are equally distributed
        size_t Mrefpre = Mrefglob * maxChunkTgtVars/Mglob;
        // ...increased by 20% to reduce the chance for required reallocation of the vectors below
        Mrefpre += Mrefpre/5;

        if (currChunk == 0) {
            // initialize chunk region with complete region first (without flanks)
            // this will be updated while reading the chunk data
            startChunkBp = startRegionBp-1; // 1-based to 0-based!
            endChunkBp = endRegionBp-1; // 1-based to 0-based!

            // prepare metadata for first chunk
            alleleFreqsCommon.reset(maxChunkTgtVars);
            if (!skipPhasing)
                cMs.reset(maxChunkTgtVars);
            // prepare target data
            for (auto &t : targets)
                t.reserve(maxChunkTgtVars);
            if (skipPhasing) {
                for (auto &tp : tgtinphase)
                    tp.reset(maxChunkTgtVars);
            }
            // prepare common reference data
            size_t capacity = roundToMultiple<size_t>(maxChunkTgtVars, UNITWORDS*sizeof(BooleanVector::data_type)*8)/8; // space for maxChunkTgtVars SNPs
            refdata = (BooleanVector::data_type*) MyMalloc::calloc(Nrefhapsmax*capacity, 1, string("refdata_c")+to_string(currChunk));
            size_t capacityT = roundToMultiple<size_t>(Nrefhapsmax, UNITWORDS*sizeof(BooleanVector::data_type)*8)/8;
            refdataT = (BooleanVector::data_type*) MyMalloc::calloc(maxChunkTgtVars*capacityT, 1, string("refdataT_c")+to_string(currChunk));
            if (!refdata || !refdataT) {
                StatusFile::addError("Not enough memory for phasing reference!");
                exit(EXIT_FAILURE);
            }
            reference.resize(Nrefhapsmax);
            auto curr_data = refdata;
            for (auto &ref : reference) {
                ref.setData(curr_data, capacity, 0);
                curr_data += capacity / sizeof(BooleanVector::data_type);
            }
            referenceT.reset(maxChunkTgtVars, maxChunkTgtVars); // also set initial size here
            auto curr_dataT = refdataT;
            for (auto &refT : referenceT) {
                refT.setData(curr_dataT, capacityT, 0);
                curr_dataT += capacityT / sizeof(BooleanVector::data_type);
            }

            currChunkOffset = 0;
            if (doImputation) // only required for imputation
                referenceFullT.reserve(Mrefpre); // don't need to prepare data memory anymore, will be done on-the-fly
            // prepare metadata for first and already for next chunk
            indexToRefFull.reset(maxChunkTgtVars);
            isImputed.reset(Mrefpre);
            nextPidx.reset(Mrefpre);
            ptgtIdx.reset(Mrefpre);

            isPhased.reset(maxChunkTgtVars);
            bcf_pout.reset(maxChunkTgtVars);

            currM = 0;
            currMref = 0;

            MLeftOv = 0;
            MrefLeftOv = 0;

            lstats.M = 0;
            lstats.Mref = 0;

        } else { // currChunk > 0

            // chunk region:
            // new start is the former end (coordinates kept equal as the chunk end could be at a multi-allelic)
            startChunkBp = endChunkBp;
            endChunkBp = endRegionBp-1; // set to region end first (1-based to 0-based), will be updated, if there are more chunks

            // DEBUG
            cerr << "LENGTHS: " << endl;
            cerr << isPhased.size() << "/"
                 << alleleFreqsCommon.size() << "/"
                 << bcf_pout.size() << "/"
                 << cMs.size() << endl;
            cerr << "tgt0: " << targets[0].size() << "/" << tgtmiss[0].size() << endl;
            for (size_t tidx = 1; tidx < Ntarget; tidx++) {
                if (targets[tidx].size() != targets[0].size())
                    cerr << "target sizes not equal!" << endl;
            }
            cerr << "MISSINGS (reverse):" << endl;
            size_t ti = 0;
            cerr << ti << ": tm  :";
            for (int i = tgtmiss[ti].size()-1; i >= 0; i--)
                cerr << "\t" << tgtmiss[ti][i];
            cerr << endl;
            // __DEBUG

            // NOTE: currentMOverlap is the number of tgt sites to take over from the last chunk as overlap

            // take over target metadata for this chunk and prepare for next chunk
            // isPhased and bcf_pout potentially carry information for non-phased tgt varianst, so they are taken care of below
            alleleFreqsCommon.keepLast(currMOverlap);
            if (!skipPhasing)
                cMs.keepLast(currMOverlap);
            // target data
            for (size_t tidx = 0; tidx < Ntarget; tidx++) {
                targets[tidx].keepLast(currMOverlap);

                // keep only missings from overlap, reduce index accordingly
                if (tgtmiss[tidx].size()) { // only if not empty
                    size_t red = currM - currMOverlap; // index reduction
                    size_t keep = 0;
                    for (int i = tgtmiss[tidx].size()-1; i >= 0; i--) { // backwards
                        if (tgtmiss[tidx][i] < red) // this missing one is not in the overlap
                            break; // further ones will neither be in the overlap
                        // else: in overlap
                        keep++;
                    }
                    tgtmiss[tidx].keepLast(keep);
                    tgtmiss[tidx].sub(red);
                }

                if (skipPhasing) {
                    tgtinphase[tidx].keepLast(currMOverlap);
                }
            }

            // DEBUG
            cerr << "MISSINGS new (reverse):" << endl;
            ti = 0;
            cerr << ti << ": tm  :";
            for (int i = tgtmiss[ti].size()-1; i >= 0; i--)
                cerr << "\t" << tgtmiss[ti][i];
            cerr << endl;
            // __DEBUG

            // reference data
            for (auto& ref : reference)
                ref.keepLast(currMOverlap);

            referenceT.keepLast(currMOverlap); // reduces the size but does not destroy the data!
            referenceT.setSize(maxChunkTgtVars); // restores the previous size and also the underlying BooleanVectors
            for (size_t i = currMOverlap; i < maxChunkTgtVars; i++) // need to clear the underlying BooleanVectors, but keep their memory assignments!
                referenceT[i].clear();

            // number of reference variants the current chunk is reduced of to keep only variants in the overlap,
            // equals former isImputed.size()-isImputedOverlap.size()
            size_t redref = indexToRefFull[indexToRefFull.size()-currMOverlap-1]+1;

            // beginning of current chunk relative to region (reference-based index)
            // -> must be the element where the first element in the overlap region points to
            currChunkOffset += redref;

            // DEBUG
            cerr << "REFLENGTHS: " << endl;
            cerr << isImputed.size() << "/"
                 << indexToRefFull.size() << "/"
                 << nextPidx.size() << "/"
                 << ptgtIdx.size() << endl;
            cerr << "refFullT.size(): " << referenceFullT.size() << " minus redref: " << (referenceFullT.size()-redref) << endl;
            // __DEBUG

            // like referenceFullT.keepLast(referenceFullT.size()-redref)
            // NOTE: referenceFullT.size() might be different to (larger than) isImputed.size() due to a potential variant swap for multi-allelics during Qref loading
            if (doImputation) { // only required for imputation
                // free data from unrequired reference haps
                for (size_t mref = 0; mref < redref; mref++) {
                    MyMalloc::free(referenceFullT[mref].getData());
                }
                // copy overlap region to beginning
                for (size_t mref = redref; mref < referenceFullT.size(); mref++) {
                    referenceFullT[mref-redref] = move(referenceFullT[mref]);
                }
                referenceFullT.resize(referenceFullT.size()-redref);
            }

            // take over reference metadata for this chunk and prepare for next chunk
            size_t keeplastref = isImputed.size()-redref; // == currMrefOverlap
            indexToRefFull.keepLast(currMOverlap);
            isImputed.keepLast(keeplastref);
            nextPidx.keepLast(keeplastref);
            ptgtIdx.keepLast(keeplastref);

            size_t keeplastisphased = isPhased.size()-nextPidx[0]; // contains overlap size + unphased sites in overlap
            isPhased.keepLast(keeplastisphased); // if outputUnphased == true -> currMOverlap + number of unphased sites in overlap
            bcf_pout.keepLast(keeplastisphased); // if outputUnphased == true -> currMOverlap + number of unphased sites in overlap

            // need to correct some index mappings now
            // DEBUG
            cerr << "idx0: " << indexToRefFull[0] << " --- "
                             << nextPidx[0] << " --- "
                             << ptgtIdx[0] << endl;
            cerr << "idxbck: " << indexToRefFull.back() << " --- "
                               << nextPidx.back() << " --- "
                               << ptgtIdx.back() << endl;
            cerr << "currMOverlap: " << currMOverlap
                 << " keeplastref: " << keeplastref
                 << " redref: " << redref
                 << " redtgt1: " << nextPidx[0]
                 << " redtgt2: " << ptgtIdx[0] << endl;
            // __DEBUG

            indexToRefFull.sub(redref);
            nextPidx.sub(nextPidx[0]);
            ptgtIdx.sub(currM-currMOverlap);

            // take over the values from the overlap from last chunk
            currM = currMOverlap; // 0 for first chunk, considered overlap for subsequent chunks
            currMref = keeplastref; // 0 for first chunk or if no imputation is done, considered overlap for subsequent chunks

            // number of vars in the left overlap:
            // the right part of the complete overlap was the part considered in last chunk.
            // now, the left part has to be considered in this chunk
            MLeftOv = currMOverlap - MRightOv;
            MrefLeftOv = keeplastref - MrefRightOv;

            // lstats should contain the stats for the current chunk without the overlap from the previous chunk,
            // so, we store the overlap here temporarily to remove it later to set the correct number
            lstats.M = currMOverlap;
            lstats.Mref = keeplastref;

        } // END currChunk > 0

        currMOverlap = 2*chunkflanksize; // set initial intended overlap size to next chunk, may be subject to change!

    } // END if(!createQRef)


    // DEBUG
    MyMalloc::printSummary(string("intermediate before reading chunk ")+to_string(currChunk+1));
//        ofstream ofs(args.outPrefix + ".memmap_c" + to_string(chunk));
//        MyMalloc::dumpMemMap(ofs);
//        ofs.close();
    // __DEBUG

    Stopwatch swrdvcf("Read data");
    stringstream ss;
    if (createQRef) {
        cout << "Creating Qref: 0%" << flush;
        ss << "Creating Qref";
    } else {
        cout << "Reading data: 0%" << flush;
        ss << "Reading data (Chunk " << (currChunk+1) << "/" << nChunks << ")";
    }
//    StatusFile::updateStatus(createQRef ? 0 : tgtlinesread/(float)Mglob, ss.str());
    StatusFile::updateStatus(0, ss.str());
    int pgb = 0; // for progress bar

    int mref_gt = 0; void *ref_gt = NULL; // will be allocated once in bcf_get_genotypes() and then reused for each marker (need void* because of htslib)
    int mtgt_gt = 0; void *tgt_gt = NULL; // will be allocated once in bcf_get_genotypes() and then reused for each marker (need void* because of htslib)

    // quick reference index set to where the last chunk stopped
    // NOTE: if the last chunk stopped at a multi-allelic and the exact match was swapped before the other alleles,
    //       qridx was behind qrefcurridxreg. This is considered here.
    size_t qridx = qrefcurridxreg - (referenceFullT.size()-isImputed.size());

    // read data SNP-wise in positional sorted order from target and reference
    // the function also takes care that we only read the requested region
//    while ((currChunk == nChunks-1 || tgtlinesread < endChunkFlankIdx[currChunk]) && bcf_sr_next_line(sr)) { // need to stay in the loop after last tgt variant if we read a VCF reference
    while ((createQRef || tgtlinesread < endChunkFlankIdx[currChunk]) && bcf_sr_next_line(sr)) {

        bcf1_t *ref = NULL;
        bcf1_t *tgt = NULL;
//        if (!loadQuickRef) {
        if (createQRef) {
            ref = bcf_sr_get_line(sr, 0); // read one line of reference, if available at current position (otherwise NULL)
            if (ref) {
//                if (createQRef && reflinesread % 16384 == 0) {
                if (reflinesread % 16384 == 0) {
                    StatusFile::updateStatus(reflinesread/(float)Mrefglob); // this is ok as long as Qrefs are created in one chunk
                    if (pgb == 3) { // print % after three dots
                        cout << (100*reflinesread/Mrefglob) << "%" << flush;
                        pgb = 0;
                    } else {
                        cout << "." << flush;
                        pgb++;
                    }
                }
                reflinesread++;
                // set the chromosome according to first reference line
                if (!setChrom) {
                    setChrom = true;
                    string chromname(bcf_hdr_id2name(ref_hdr, ref->rid)); // get the chromosome name, as written in the reference file
                    string chromlittmp;
                    int chromtmp = chrNameToNumber(chromname, chromlittmp, &chrliterally); // get the chromosome number
                    if ((chrom > 0 && chrom != chromtmp) || (chrom == 0 && chromlit.compare(chromlittmp))) {
                        cout << endl;
                        StatusFile::addError("Chromosome mismatch between reference filename and content.");
                        exit(EXIT_FAILURE);
                    }
                    chromlit = chromlittmp; // just to be sure to take the notation from inside the reference file (the only allowed difference is the preceding "chr")
                    if (currChunk == 0) {
                        initInfoFile();
                        if (doImputation)
                            initImpHeader();
                    }
                }
            }
        }
        if (!createQRef) {
            tgt = bcf_sr_get_line(sr, loadQuickRef ? 0 : 1); // read one line of target, if available at current position (otherwise NULL)
            if (tgt) {
                if (tgtlinesread % 1024 == 0) { // TODO as the chunk size is (approximately) known, we can use this to better display the progress (i.e. change the modulus)
                    float progress = currChunk == 0 ? (tgtlinesread/(float)endChunkFlankIdx[0]) : ((tgtlinesread - endChunkFlankIdx[currChunk-1])/(float)(endChunkFlankIdx[currChunk]-endChunkFlankIdx[currChunk-1]));
                    int progresspercent = currChunk == 0 ? (100*tgtlinesread/endChunkFlankIdx[0]) : (100*(tgtlinesread - endChunkFlankIdx[currChunk-1])/(endChunkFlankIdx[currChunk]-endChunkFlankIdx[currChunk-1]));
                    StatusFile::updateStatus(progress);
                    if (pgb == 3) { // print % after three dots
                        cout << progresspercent << "%" << flush;
                        pgb = 0;
                    } else {
                        cout << "." << flush;
                        pgb++;
                    }
                }

                // check, how we get along with our memory:
                // if we already used up ~97% of our available memory, we stop this chunk now
                // NOTE: we still load the reference variants between the last and this tgt site
                // DEBUG
                if (MyMalloc::getCurrentAlloced() > maxchunkmem-maxchunkmem/32) {
//                if (MyMalloc::getCurrentAlloced() > maxchunkmem-maxchunkmem/4) { // see if increasing nchunks works

                    // TODO maybe adjust currMOverlap here, if the chunk is too small
                    // -> but this needs to be considered for the calculation of the total number of chunks later!!
                    // -> and for the calculation of MRightOv/MrefRightOv!!
                    int64_t reduction = endChunkFlankIdx[currChunk] - tgtlinesread - 1; // difference to the previously calculated end of chunk
                    startChunkFlankIdx[currChunk+1] -= reduction;
                    startChunkIdx[currChunk+1] -= reduction;
                    endChunkFlankIdx[currChunk] -= reduction; // this will also stop reading after this line

                    // DEBUG
                    if (reduction)
                        cerr << "--- REDUCED chunk by " << reduction << " sites." << endl;
                }

                tgtlinesread++;

                // check, if the target contains the same chromosome as the chosen reference
                if (!checkedChrom) {
                    checkedChrom = true;
                    string chromname(bcf_hdr_id2name(tgt_hdr, tgt->rid)); // get the chromosome name, as written in the target file
                    string chromlittmp;
                    bool chrliterally_tmp;
                    int chromtmp = chrNameToNumber(chromname, chromlittmp, &chrliterally_tmp); // get the chromosome number
                    if ((chrom > 0 && chrom != chromtmp) || (chrom == 0 && chromlit.compare(chromlittmp))) {
                        cout << endl;
                        StatusFile::addError("Chromosome mismatch between reference and target.");
                        exit(EXIT_FAILURE);
                    }
                    // we need to init the imputation header here only if we load a Qref or if the first line is a tgt line
                    if (!setChrom || loadQuickRef) {
                        setChrom = true; // if not loading a Qref, we take the style of the target file only in the rare case if the first line is from the target
                        chromlit = chromlittmp; // just to be sure to take the notation from inside the target file
                        chrliterally = chrliterally_tmp;
                        if (currChunk == 0) {
                            initInfoFile();
                            if (doImputation)
                                initImpHeader();
                        }
                    }
                }
            } // END if (tgt)
        } // END if (!createQref)

        if (useExclude && bcf_sr_get_line(sr, loadQuickRef ? 1 : 2)) { // SNP is in exclude file -> skip
          if (tgt) {
              lstats.Mexclude++;
              if (outputUnphased) { // this SNP is explicitly excluded so it will be untouched for phasing, but if the user decides to "outputUnphased" it will be copied to the output files
                  bcf_pout.push_back(bcf_dup(tgt));
                  isPhased.push_back(false);
                  lstats.Munphased++;
              }
              addToInfoFileExcluded(tgt, "exclude file");
          }
          continue;
        }

        // if we read a quick reference, tgt is set here and needs to be processed.
        // so, we forward our qref cursor either to a common variant or to the next variant (in which case the current variant is target-only).
        size_t qfoundidx = 0, qrold = qridx;
        bool qrswapped = false, qrflipped = false, qrrefalterror = false;
        bool qfound = loadQuickRef ? loadAndFindInQref(tgt, qridx, qrswapped, qrflipped, qrrefalterror) : false;
        if (loadQuickRef) {
            if (qfound) { // if we found the tgt variant, and there were no errors, we continue with the next one next time
                qfoundidx = qridx;
                if (!qrrefalterror)
                    qridx++;
            }
            // anyway, we need to process the variants now
            for (size_t qidx = qrold; qidx < qridx; qidx++) {
                isImputed.push_back(true); // this is the default and may change below if we found a common SNP here that will be phased
                nextPidx.push_back(bcf_pout.size());
                ptgtIdx.push_back(currM);
            }
            currMref += qridx - qrold;
        }

        if ((!ref && !qfound) || (tgt && !loadQuickRef && ref->n_allele == 1)) { // SNP is not in reference (thus, it has to be in target only), or ref is mono but SNP is also in target
//            if (tgt->n_allele > 1) // report if polymorphic in target -> should be counted even if monomorphic
            if (ref) { // means ref->n_allele == 1
                lstats.MmonomorphicRefTgt++;
                addToInfoFileExcluded(tgt, "monomorphic reference");
            } else {
                lstats.MtargetOnly++;
                addToInfoFileExcluded(tgt, "target only");
            }
            if (outputUnphased) { // this SNP cannot be used for phasing and imputation, but if the user decides to "outputUnphased" it will be copied to the phased output file
                bcf_pout.push_back(bcf_dup(tgt));
                isPhased.push_back(false);
                lstats.Munphased++;
            }
            continue;
        }

        if (!tgt) { // SNP is in reference, but not in target -> need to keep it for imputation (if desired) (usual case if we create a qref, but not the case if we load a qref)

            lstats.MrefOnly++;
            if (doImputation || createQRef) {
                int numMissing, numUnphased;
                if (ref->n_allele <= 2) {

                    processReferenceOnlySNP(ref, &ref_gt, &mref_gt, numMissing, numUnphased, false);
                    if (numMissing)
                        lstats.MmissingRefOnly++;
                    if (numUnphased)
                        lstats.MunphasedRefOnly++;
                    lstats.GmissingRefOnly += numMissing;
                    lstats.GunphasedRefOnly += numUnphased;

                } else {
                    lstats.MmultiAllelicRefOnly++;

                    if (!excludeMultiAllRef || createQRef) { // multi-allelic: split and keep only if the user does not exclude multi-allelic SNPs
                        vector<bcf1_t*> masplits;
                        splitSNP(ref, &ref_gt, &mref_gt, masplits, lstats.MmultiSplittedRefOnly, lstats.GmultiFilledRefOnly);
                        processReferenceOnlySNP(ref, &ref_gt, &mref_gt, numMissing, numUnphased, true);
                        if (numMissing)
                            lstats.MmissingRefOnly++;
                        if (numUnphased)
                            lstats.MunphasedRefOnly++;
                        lstats.GmissingRefOnly += numMissing;
                        lstats.GunphasedRefOnly += numUnphased;
                        // store remaining splits
                        for (auto splitit = masplits.begin(); splitit != masplits.end(); splitit++) {
                            if (doImputation || createQRef) {
                                lstats.MrefOnly++;
                                processReferenceOnlySNP(*splitit, &ref_gt, &mref_gt, numMissing, numUnphased, true);
                                if (numMissing)
                                    lstats.MmissingRefOnly++;
                                if (numUnphased)
                                    lstats.MunphasedRefOnly++;
                                lstats.GmissingRefOnly += numMissing;
                                lstats.GunphasedRefOnly += numUnphased;
                            }
                        }
                        // need to clean up the generated splits
                        for (bcf1_t* s : masplits)
                            bcf_destroy(s);
                    }
                }
            }
            continue;
        }

        // Here, we have a common target and ref!

        // ***
        // if we create a qref file, we will not reach this point!
        // ***

        // deal with multi-allelic and monomorphic markers
        // drop multi-allelic target markers
        if (tgt->n_allele > 2) {
            lstats.MmultiAllelicTgt++;
            if (outputUnphased) {
                bcf_pout.push_back(bcf_dup(tgt));
                isPhased.push_back(false);
                lstats.Munphased++;
            }
            addToInfoFileExcluded(tgt, "multi-allelic target");
            // this SNP is excluded in imputation as well, so no output
            // (but will be imputed if loaded a qref...)
            continue;
        }

        // drop multi-allelic reference markers if desired by user
        if (excludeMultiAllRef) {
//            if ((!loadQuickRef && ref->n_allele > 2)
//                    || (loadQuickRef && multiAllFlagsFullRefRegion[qfoundidx])) {
            if (multiAllFlagsFullRefRegion[qfoundidx]) {
                lstats.MmultiAllelicRefTgt++;
                if (outputUnphased) {
                    bcf_pout.push_back(bcf_dup(tgt));
                    isPhased.push_back(false);
                    lstats.Munphased++;
                }
                addToInfoFileExcluded(tgt, "multi-allelic reference");
                // this SNP is excluded in imputation as well, so no output
                continue;
            }
        }

//        // if not done so far, we need to unpack the records here:
//        // for accessing ref->d.allele we need to unpack up to ALT at least (BCF_UN_STR),
//        // (in order to access the genotypes, it is called again later with BCF_UN_FMT
//        // within bcf_get_genotypes())
//        if (!loadQuickRef)
//            bcf_unpack(ref, BCF_UN_STR);

        // we can fetch the genotypes already here, which also unpacks the d.allele fields
        int ntgt_gt = bcf_get_genotypes(tgt_hdr, tgt, &tgt_gt, &mtgt_gt); // calls bcf_unpack() within

        // for info file: save original alleles
        string tref(tgt->d.allele[0]);
        string talt;
        if (tgt->n_allele > 1)
            talt = string(tgt->d.allele[1]);

        // preserve monomorphic markers if not monomorphic in the reference panel
        if (tgt->n_allele < 2) {
//            if (loadQuickRef) {
            string allref = qrswapped ? allelesFullRefRegion[2*qfoundidx+1] : allelesFullRefRegion[2*qfoundidx];
            string allalt = qrswapped ? allelesFullRefRegion[2*qfoundidx] : allelesFullRefRegion[2*qfoundidx+1];
            if (qrflipped) {
                allref = reverseComplement(allref.c_str());
                allalt = reverseComplement(allalt.c_str());
            }
            bcf_update_alleles_str(tgt_hdr, tgt, allref.append(",").append(allalt).c_str());
//            } else {
//                // check which allele is ours
//                unsigned int a = 0;
//                bool strandflip = false;
//                bool found = true; // we expect to find the allele
//                while (a < ref->n_allele && strcmp(ref->d.allele[a], tgt->d.allele[0]) != 0)
//                    a++;
//                if (a == ref->n_allele) { // not found: strand flip?
//                    a = 0;
//                    while (a < ref->n_allele && reverseComplement(ref->d.allele[a]).compare(tgt->d.allele[0]))
//                        a++;
//                    if (a == ref->n_allele) // not found
//                        found = false;
//                    else // strand flip
//                        strandflip = true;
//                }
//                if (found) {
//                    if (!strandflip) {
//                        if (a) // swap alleles if the common one was not the reference
//                            bcf_update_alleles_str(tgt_hdr, tgt, string(ref->d.allele[a]).append(".").append(ref->d.allele[0]).c_str());
//                        else
//                            bcf_update_alleles_str(tgt_hdr, tgt, string(ref->d.allele[0]).append(".").append(ref->d.allele[1]).c_str());
//                    } else { // strand flip
//                        if (a)
//                            bcf_update_alleles_str(tgt_hdr, tgt, string(reverseComplement(ref->d.allele[a])).append(".").append(reverseComplement(ref->d.allele[0])).c_str());
//                        else
//                            bcf_update_alleles_str(tgt_hdr, tgt, string(reverseComplement(ref->d.allele[0])).append(".").append(reverseComplement(ref->d.allele[1])).c_str());
//                    }
//                } else
//                    bcf_update_alleles_str(tgt_hdr, tgt, string(tgt->d.allele[0]).append(",.").c_str()); // make bi-allelic with one allele missing. will later result in refalterror anyway
//            }
        }

        // check for REF/ALT errors
//        bool refaltswap = loadQuickRef ? qrswapped : false;
//        bool strandflip = loadQuickRef ? qrflipped : false;
//        bool refalterror = loadQuickRef ? (qrrefalterror || (qrswapped && !allowRefAltSwap) || (qrflipped && !allowStrandFlip)) : true; // default to true, if we find a match, we set it to false
        bool refaltswap = qrswapped;
        bool strandflip = qrflipped;
        bool refalterror = qrrefalterror || (qrswapped && !allowRefAltSwap) || (qrflipped && !allowStrandFlip);
//        if (!loadQuickRef) {
//            // test through all possibilities in order, stop if we find a match:
//            // 1. direct match
//            if (strcmp(tgt->d.allele[0], ref->d.allele[0]) == 0) { // equal reference alleles
//                if (strcmp(tgt->d.allele[1], ref->d.allele[1]) == 0) { // equal first alternative alleles
//                    refalterror = false;
//                } else if (ref->n_allele > 2) { // if ref is multi-allelic, need to find the alt allele
//                    for (int i = 2; i < ref->n_allele; i++) {
//                        if (strcmp(tgt->d.allele[1], ref->d.allele[i]) == 0) {
//                            refalterror = false;
//                            break;
//                        }
//                    }
//                }
//            } // NOTE: we will stop here if we've already found an equal SNP here as refalterror was already set to false then
//            // alternatives are only possible if tgt is a SNP
//            if (bcf_is_snp(tgt)) {
//                // 2. ref/alt swapped
//                if (refalterror && allowRefAltSwap && strcmp(tgt->d.allele[1], ref->d.allele[0]) == 0) { // equal alternative target allele and ref reference allele
//                    if (strcmp(tgt->d.allele[0], ref->d.allele[1]) == 0) { // target ref matches first alternative reference allele
//                        refalterror = false;
//                        refaltswap = true;
//                    } else if (ref->n_allele > 2) { // if ref is multi-allelic, need to find the alt allele
//                        for (int i = 2; i < ref->n_allele; i++) {
//                            if (strcmp(tgt->d.allele[0], ref->d.allele[i]) == 0) {
//                                refalterror = false;
//                                refaltswap = true;
//                                break;
//                            }
//                        }
//                    }
//                }
//                // 3. strand flipped
//                if (refalterror && allowStrandFlip && reverseComplement(tgt->d.allele[0]).compare(ref->d.allele[0]) == 0) { // equal reference alleles if reverse complemented
//                    if (reverseComplement(tgt->d.allele[1]).compare(ref->d.allele[1]) == 0) { // equal first alternative alleles if reverse complemented
//                        refalterror = false;
//                        strandflip = true;
//                    } else if (ref->n_allele > 2) { // if ref is multi-allelic, need to find the alt allele
//                        for (int i = 2; i < ref->n_allele; i++) {
//                            if (reverseComplement(tgt->d.allele[1]).compare(ref->d.allele[i]) == 0) {
//                                refalterror = false;
//                                strandflip = true;
//                                break;
//                            }
//                        }
//                    }
//                }
//                // 4. ref/alt swapped and strand flipped
//                if (refalterror && allowRefAltSwap && allowStrandFlip && reverseComplement(tgt->d.allele[1]).compare(ref->d.allele[0]) == 0) { // equal reference alleles if reverse complemented
//                    if (reverseComplement(tgt->d.allele[0]).compare(ref->d.allele[1]) == 0) { // target ref matches first alternative reference allele if reverse complemented
//                        refalterror = false;
//                        refaltswap = true;
//                        strandflip = true;
//                    } else if (ref->n_allele > 2) { // if ref is multi-allelic, need to find the alt allele
//                        for (int i = 2; i < ref->n_allele; i++) {
//                            if (reverseComplement(tgt->d.allele[0]).compare(ref->d.allele[i]) == 0) {
//                                refalterror = false;
//                                refaltswap = true;
//                                strandflip = true;
//                                break;
//                            }
//                        }
//                    }
//                }
//            } // END if is SNP?
//        } // END !loadQref

        // skip on ref/alt error
        if (refalterror) {
            lstats.MrefAltError++;

            if (outputUnphased) {
                bcf_pout.push_back(bcf_dup(tgt));
                isPhased.push_back(false);
                lstats.Munphased++;
            }
            addToInfoFileExcluded(tgt->pos, tgt->d.id, tref, talt, "ref/alt error");
            // this SNP is excluded in imputation as well, so no output
            // (but will be imputed if loaded a qref...)
            continue;
        }

        if (refaltswap && strandflip) {
            // replace allele strings by their reverse complements and swap
            bcf_update_alleles_str(tgt_hdr, tgt, string(reverseComplement(tgt->d.allele[1]).append(",").append(reverseComplement(tgt->d.allele[0]))).c_str());
            lstats.numRefAltSwapAndStrandFlip++;
        } else if (refaltswap) {
            // swap allele strings, genotypes are interpreted correspondingly below when processing target SNP
            swap(tgt->d.allele[0], tgt->d.allele[1]);
            bcf_update_alleles(tgt_hdr, tgt, (const char**) (tgt->d.allele), 2);
            lstats.numRefAltSwaps++;
        } else if (strandflip) {
            // replace allele strings by their reverse complements and swap
            bcf_update_alleles_str(tgt_hdr, tgt, string(reverseComplement(tgt->d.allele[0]).append(",").append(reverseComplement(tgt->d.allele[1]))).c_str());
            lstats.numStrandFlips++;
        }

        // SNP passes checks

        // append base pair coordinate to chrBps
        chrBpsReg.push_back((uint64_t)(tgt->pos + 1));
        if (!skipPhasing) { // only required (and available) for phasing
            cMs.push_back(mapint.interp(chrBpsReg.back()));
        }

//        // split reference SNP into several SNPs if multi-allelic (if we wanted to exclude multi-allelics, we would have done already)
//        vector<bcf1_t*> masplits;
//        if (!loadQuickRef) {
//            if (ref->n_allele > 2) {
//                lstats.MmultiAllelicRefTgt++;
//                int numMissing, numUnphased;
//                splitSNP(ref, &ref_gt, &mref_gt, masplits, lstats.MmultiSplittedTgt, lstats.GmultiFilledTgt);
//
//                auto splitit = masplits.begin();
//                while( splitit != masplits.end() && strcmp(tgt->d.allele[1], ref->d.allele[1]) != 0 ) {
//                    // current ref is not the one we need for phasing here, so store only for imputation and try next
//                    if (doImputation) {
//                        lstats.MrefOnly++;
//                        processReferenceOnlySNP(ref, &ref_gt, &mref_gt, numMissing, numUnphased, inOverlap, true);
//                        if (numMissing)
//                            lstats.MmissingRefOnly++;
//                        if (numUnphased)
//                            lstats.MunphasedRefOnly++;
//                        lstats.GmissingRefOnly += numMissing;
//                        lstats.GunphasedRefOnly += numUnphased;
//                    }
//                    ref = *splitit;
//                    bcf_unpack(ref, BCF_UN_STR);
//                    splitit++;
//                }
//                // store remaining splits (will be listed before the one used for imputation/phasing in the output)
//                for (; splitit != masplits.end(); splitit++) {
//                    if (doImputation) {
//                        lstats.MrefOnly++;
//                        processReferenceOnlySNP(*splitit, &ref_gt, &mref_gt, numMissing, numUnphased, inOverlap, true);
//                        if (numMissing)
//                            lstats.MmissingRefOnly++;
//                        if (numUnphased)
//                            lstats.MunphasedRefOnly++;
//                        lstats.GmissingRefOnly += numMissing;
//                        lstats.GunphasedRefOnly += numUnphased;
//                    }
//                }
//                multiAllFlagsFullRefRegion.push_back_withPreInit(true);
//                MrefMultiAllreg++;
//            } else { // n_allele <= 2
//                multiAllFlagsFullRefRegion.push_back_withPreInit(false);
//            }
//        } else { // loadQuickRef
        if (multiAllFlagsFullRefRegion[qfoundidx])
            lstats.MmultiAllelicRefTgt++; // just for the statistics
//        }

        // include SNP in reference data for phasing and, if desired, in the full data for imputation

//        int numMissing = 0, numUnphased = 0;
        int numMissing = 0;
        float af = 0.0;
//        if (loadQuickRef) {
        // copy ref information for phasing
        af = alleleFreqsFullRefRegion[qfoundidx];
        referenceT[currM].setSize(Nrefhaps); // size will be the number of haps to add
        for (size_t i = 0; i < Nrefhaps; i++) {
            reference[i].push_back_withPreInit(referenceFullT[qfoundidx-currChunkOffset][i]);
            referenceT[currM].setWithPreInit(i, referenceFullT[qfoundidx-currChunkOffset][i]);
        }
//        } else { // !loadQuickRef
//            processReferenceSNP(Nref, ref, &ref_gt, &mref_gt, inOverlap, numMissing, numUnphased, af, true);
//
//            if (numMissing)
//                lstats.MwithMissingRef++;
//            if (numUnphased)
//                lstats.MwithUnphasedRef++;
//            lstats.GmissingRef += numMissing;
//            lstats.GunphasedRef += numUnphased;
//        }

        // process target genotypes: append Ntarget entries (0/1/2/9) to genosTarget[]
        processTargetSNP(Ntarget, ntgt_gt, reinterpret_cast<int*>(tgt_gt), refaltswap, numMissing, tgt->pos);
        lstats.GmissingTarget += numMissing;

//        // keep the record's information
//        if (doImputation && !loadQuickRef) {
//            alleleFreqsFullRefRegion.push_back(af);
//            positionsFullRefRegion.push_back(ref->pos);
//            allelesFullRefRegion.emplace_back(ref->d.allele[0]);
//            allelesFullRefRegion.emplace_back(ref->d.allele[1]);
//            variantIDsFullRefRegion.emplace_back(ref->d.id);
//
//            isImputed.push_back(false);
//            indexToRefFull.push_back(currMref);
//            nextPidx.push_back(bcf_pout.size()); // common site, index will point to element being inserted now
//            ptgtIdx.push_back(currM);
//            currMref++;
//            if (inOverlap) {
//                isImputedOverlap.push_back(false);
//                indexToRefFullOverlap.push_back(currMrefOverlap);
//                nextPidxOverlap.push_back(bcf_poutOverlap.size()); // common site, index will point to element being inserted now
//                ptgtIdxOverlap.push_back(currMOverlap);
//                currMrefOverlap++;
//            }
//        }
//        if (loadQuickRef) { // loaded quick ref
        // we need to change the imputation flag since it was initialized with "true"
        isImputed.set_back(false);
        indexToRefFull.push_back(isImputed.size()-1);
        // nextPidx and ptgtIdx was already taken care of after loading the Qref
//        }

        // remove genotypes (will be phased)
        bcf_update_genotypes(tgt_hdr, tgt, NULL, 0);
        isPhased.push_back(true);
        bcf_pout.push_back(bcf_dup(tgt));
        alleleFreqsCommon.push_back(af);
        currM++;
        string infoexpl;
        if (refaltswap && strandflip)
            infoexpl = "ref/alt swap + strand flip";
        else if (strandflip)
            infoexpl = "strand flip";
        else if (refaltswap)
            infoexpl = "ref/alt swap";
//        if (loadQuickRef)
        addToInfoFileIncluded(tgt->pos, tgt->d.id, tref, talt, variantIDsFullRefRegion[qfoundidx], allelesFullRefRegion[2*qfoundidx], allelesFullRefRegion[2*qfoundidx+1], infoexpl);
//        else
//            addToInfoFileIncluded(tgt->pos, tgt->d.id, tref, talt, ref->d.id, ref->d.allele[0], ref->d.allele[1], infoexpl);

//        // need to clean up potentially generated splits -> if ref points to one of these splits, it is going to be destroyed here! don't use it afterwards!
//        for (bcf1_t* s : masplits)
//            bcf_destroy(s);

    } // END while (reading chunk)

    // allocated by HTSlib
    free(ref_gt);
    free(tgt_gt);

    if (loadQuickRef) {
        if (tgtlinesread == Mglob) { // if we are at the end of the last chunk, the last ref-only variants have to be processed
            // no overlap
            MRightOv = 0;
            MrefRightOv = 0;

            // read remaining reference vars
            while (qrefcurridxreg < Mrefreg) {
                // TODO we have to check the memory here and break up if there are too many reference variants left!
                qRefLoadNextVariant();
                isImputed.push_back(true);
                nextPidx.push_back(bcf_pout.size());
                ptgtIdx.push_back(currM);
                currMref++;
            }

            // set end of chunk to last reference var (inclusive)
            endChunkBp = positionsFullRefRegion.back();
            // DEBUG
            cerr << "endChunkBp: " << endChunkBp << endl;

            // no more chunks
            startChunkFlankIdx[currChunk+1] = endChunkFlankIdx[currChunk];
            startChunkIdx[currChunk+1] = endChunkFlankIdx[currChunk];
            currMOverlap = 0;

            // set the final number of chunks
            nChunks = currChunk + 1;

        } else { // there will be another chunk
            // overlap sizes
            MRightOv = currMOverlap/2; // number of tgt variants in the overlap to next chunk (chunkflanksize for now)
            MrefRightOv = currMref - indexToRefFull[currM-MRightOv-1] - 1; // number of ref vars in overlap to next chunk

            // set end of chunk to bp position of last common var in this chunk
            endChunkBp = positionsFullRefRegion[currChunkOffset+indexToRefFull[currM-chunkflanksize-1]]; // inclusive!

            // DEBUG
            cerr << "endChunkBp: " << endChunkBp << " -- " << bcf_pout[currM-chunkflanksize-1]->pos << " / " << bcf_pout[currM-chunkflanksize]->pos << " / " << bcf_pout[currM-chunkflanksize+1]->pos << endl;

            // prepare preliminary limits
            size_t nextchunkend = min(Mglob, startChunkFlankIdx[currChunk+1] + maxChunkTgtVars);
            endChunkFlankIdx.push_back(nextchunkend);
            startChunkIdx.push_back(nextchunkend - chunkflanksize);
            startChunkFlankIdx.push_back(nextchunkend - 2*chunkflanksize);

            // depending on the potentially reduced chunk size, we can correct the number of chunks now
            int remChunks = nChunks-currChunk-1; // how many chunks are still planned after this chunk
            int64_t remvars = Mglob - startChunkFlankIdx[currChunk+1]; // variants that still have to be processed (including the current overlap)
            remvars -= remChunks * maxChunkTgtVars; // how many vars are spanned by this chunk...
            if (remChunks) // ...with consideration of overlaps
                remvars += (nChunks-currChunk-2)*2*chunkflanksize;
            if (remvars > 0)
                nChunks++;
            // DEBUG
            cerr << "Remvars: " << remvars << endl;
            // _DEBUG
        }
    }

    // DEBUG
    cerr << "qrefcurridxglob: " << qrefcurridxglob << " qrefcurridxreg: " << qrefcurridxreg << endl;

    // set M and Mref for the current chunk
    M = currM;
    Mref = currMref;

    lstats.M = currM - lstats.M; // we have to remove the overlap size from previous chunk as it is contained in currM
    lstats.Mref = currMref - lstats.Mref; // we have to remove the overlap size from previous chunk as it is contained in currMref
    lstats.MrefMultiAllreg = MrefMultiAllreg;

    // reduce transposed reference according to the number of variants we inserted
    referenceT.setSize(M);

    // set the map for haploid samples
    if (currChunk == 0 && !createQRef) {
        // reference
//        if (!loadQuickRef) { // otherwise already set
////            if (chrom == CHRX || chrom == CHRY) { // go through the flags if we are on X or Y
//                for (size_t i = 0; i < Nrefhapsmax/2; i++) {
//                    haploidsRefMap.push_back(2*i);
//                    if (!haploidsRef[i]) // diploid
//                        haploidsRefMap.push_back(2*i+1);
//                }
////            } else { // diploid chromosome: identity
////                for (size_t i = 0; i < Nrefhapsmax; i++)
////                    haploidsRefMap.push_back(i);
////            }
//        }
        // target
//        if (chrom == CHRX || chrom == CHRY) { // go through the flags if we are on X or Y
            for (size_t i = 0; i < Ntarget; i++) {
                haploidsTgtMap.push_back(2*i);
                if (!haploidsTgt[i]) // diploid
                    haploidsTgtMap.push_back(2*i+1);
            }
//        } else { // diploid chromosome: identity
//            for (size_t i = 0; i < 2*Ntarget; i++)
//                haploidsTgtMap.push_back(i);
//        }
    }

    // convert variant IDs to chr:pos:ref:alt if missing
    if (noMissingIDs) {
        for (size_t i = 0; i < variantIDsFullRefRegion.size(); i++) {
            // variant IDs starting with '.' are considered missing
            if (variantIDsFullRefRegion[i][0] == '.') { // empty variant IDs are not allowed
                variantIDsFullRefRegion[i] = chromlit + ":" + to_string(positionsFullRefRegion[i]+1) + ":" + allelesFullRefRegion[2*i] + ":" + allelesFullRefRegion[2*i+1];
            }
        }
    }

    // calculate the number of output sites for each file (left and right parts in the overlaps are not considered as they won't be written here)
    num_sites_per_file[0] = MrefLeftOv; // left overlap
    num_sites_per_file.back() = MrefRightOv; // right overlap
    // number of remaining sites per file
    size_t mref_rem = Mref - MrefLeftOv - MrefRightOv;
    size_t nspf = roundToMultiple(mref_rem, (size_t) num_files) / num_files;
    for (unsigned f = 0; f < num_files; f++) {
        num_sites_per_file[f+1] = nspf;
        if (mref_rem % num_files && f >= mref_rem % num_files) {
            // evenly distribute
            num_sites_per_file[f+1]--;
        }
    }

    swrdvcf.stop();

    // DEBUG
    cerr << "sf/s/e: " << startChunkFlankIdx[currChunk] << "/" << startChunkIdx[currChunk] << "/" << endChunkFlankIdx[currChunk]
         << " size/sizef: " <<  endChunkFlankIdx[currChunk] - startChunkIdx[currChunk] << "/" << endChunkFlankIdx[currChunk] - startChunkFlankIdx[currChunk]
         << " tgtlinesread: " << tgtlinesread
         << " sfn/sn: " << startChunkFlankIdx[currChunk+1] << "/" << startChunkIdx[currChunk+1]
         << endl;
    cerr << "currChunkOffset: " << currChunkOffset << " currM: " << M << " currMref: " << Mref << endl;
    cerr << "MLeftOv: " << MLeftOv << " MRightOv: " << MRightOv << endl;
    cerr << "MrefLeftOv: " << MrefLeftOv << " MrefRightOv: " << MrefRightOv << endl;
    cerr << "startBp: " << startChunkBp << " endBp: " << endChunkBp << endl;
    cerr << "positions:" << endl;

    cerr << " L: " << positionsFullRefRegion[currChunkOffset] << "/(" // first ref variant loaded (common or not?)
                   << positionsFullRefRegion[currChunkOffset+MrefLeftOv] << "/" // first ref variant in this chunk (common or not?)
                   << positionsFullRefRegion[currChunkOffset+(MLeftOv?indexToRefFull[MLeftOv-1]:0)] << "/" // last common variant in previous chunk
                   << positionsFullRefRegion[currChunkOffset+indexToRefFull[MLeftOv]] << ")" << endl; // first common variant in this chunk

    cerr << " R: " << positionsFullRefRegion[currChunkOffset+indexToRefFull[M-currMOverlap-1]] << "/" // last common not in overlap
                   << positionsFullRefRegion[currChunkOffset+indexToRefFull[M-currMOverlap]] << "/(" // first common in overlap
                   << positionsFullRefRegion[currChunkOffset+Mref-MrefRightOv-1] << "/" // last ref in this chunk (middle of overlap)
                   << positionsFullRefRegion[currChunkOffset+Mref-MrefRightOv] << "/" // first ref in next chunk (middle of overlap)
                   << positionsFullRefRegion[currChunkOffset+indexToRefFull[M-MRightOv-1]] << "/" // last common in this chunk (middle of overlap)
                   << positionsFullRefRegion[currChunkOffset+indexToRefFull[M-MRightOv]] << ")/(" // first common in next chunk (middle of overlap)
                   << positionsFullRefRegion[currChunkOffset+indexToRefFull[M-1]] << "/" // last common loaded
                   << positionsFullRefRegion[currChunkOffset+Mref-1] << ")" << endl; // last ref loaded (common if not last chunk, else last ref at all)




    if (pgb == 0) // just printed "xx%"
        cout << ".";
    cout << "100%";
    if (currChunk+1 != nChunks)
        cout << " (chunk end)" << endl;
    cout << endl;

    if (!createQRef) { // statistics for normal run

//        if (loadQuickRef)
        lstats.MmultiAllelicRefOnly = MrefMultiAllreg - lstats.MmultiAllelicRefTgt;

        // number of reference variants the current chunk is reduced of to keep only variants in the overlap,
        // equals former isImputed.size()-isImputedOverlap.size()
        size_t redref = indexToRefFull[indexToRefFull.size()-currMOverlap-1]+1;
        size_t keeplastref = indexToRefFull.back()+1-redref; // == currMrefOverlap

        stringstream stmp;
        stmp << "<p class='pinfo'><b>" << M << " variants in both target and reference are used for phasing.</b>";
        if (currChunk+1 < nChunks)
            stmp << "<br>\n  (Of these are " << currMOverlap << " variants in the overlap to next chunk.)";
        stmp << "</p>";
        size_t Mreftmp = loadQuickRef ? Mref-MrefMultiAllreg : Mref;
        if (doImputation) {
            stmp << "\n<p class='pinfo'>" << (Mreftmp-M) << " variants exclusively in reference will be imputed.<br>\n";
            if (currChunk+1 < nChunks)
                stmp << "  (Of these are " << (keeplastref - currMOverlap) << " variants in the overlap to next chunk.)<br>\n";
            stmp << "<b>Imputation output will contain " << Mreftmp << " variants.</b>";
            if (currChunk+1 < nChunks)
                stmp << "<br>\n  (" << keeplastref << " in the overlap to next chunk.)";
             stmp << "</p>";
        }
        StatusFile::addInfo(stmp.str(), false);

        // for YAML messages, placed only if yaml is set
        stringstream yamlinfo;
        yamlinfo << "\n";
        yamlinfo << "    Index: " << (currChunk+1) << "\n";
        yamlinfo << "    Common variants: " << M << "\n";
        if (currChunk+1 < nChunks)
            yamlinfo << "    Common variants in overlap: " << currMOverlap << "\n";
        if (doImputation) {
            yamlinfo << "    Exclusive reference variants: " << (Mreftmp - M) << "\n";
            if (currChunk+1 < nChunks)
                yamlinfo << "    Exclusive reference variants in overlap: " << (keeplastref - currMOverlap) << "\n";
            yamlinfo << "    Imputation output variants: " << Mreftmp << "\n";
            if (currChunk+1 < nChunks)
                yamlinfo << "    Imputation output variants in overlap: " << keeplastref << "\n";
        }

        // determine SNP rate (only if we do phasing)
        if (!skipPhasing && nHaploidsTgt != Ntarget) {
            size_t physRange = M <= 1 ? 0 : chrBpsReg.back()-chrBpsReg[chrBpsReg.size()-cMs.size()];
            fp_type cMrange = M <= 1 ? 0 : cMs.back() - cMs[0];
            // store first cM value for summary statistics
            if (currChunk == 0)
                cM0 = cMs[0];

            StatusFile::addInfo("<table>", false);
            StatusFile::addInfo("<tr><td>Physical distance range:        </td><td>" + to_string(physRange) + " base pairs</td></tr>", false);
            StatusFile::addInfo("<tr><td>Genetic distance range:         </td><td>" + to_string(cMrange) + " cM</td></tr>", false);
            StatusFile::addInfo("<tr><td>Average #SNPs per cM in target: </td><td>" + to_string((int)(M/cMrange + 0.5)) + "</td></tr></table>", false);

            yamlinfo << "    Physical distance range: " << physRange << " bp\n";
            yamlinfo << "    Genetic distance range: " << cMrange << " cM\n";
            yamlinfo << "    Average SNPs per cM in target: " << (int)(M/cMrange + 0.5) << "\n";
        }

        StatusFile::addInfo("<p class='pinfo'>", false);
        if (lstats.numRefAltSwapAndStrandFlip) {
            StatusFile::addInfo("  REF/ALT were swapped AND strands were flipped in " + to_string(lstats.numRefAltSwapAndStrandFlip) + " variants.<br>", false);
            yamlinfo << "    REF ALT swap strand flips: " << lstats.numRefAltSwapAndStrandFlip << "\n";
        }
        if (lstats.numRefAltSwaps) {
            StatusFile::addInfo("  REF/ALT were swapped in " + to_string(lstats.numRefAltSwaps) + " variants.<br>", false);
            yamlinfo << "    REF ALT swaps: " << lstats.numRefAltSwaps << "\n";
        }
        if (lstats.numStrandFlips) {
            StatusFile::addInfo("  Strands were flipped in " + to_string(lstats.numStrandFlips) + " variants.<br>", false);
            yamlinfo << "    Strand flips: " << lstats.numStrandFlips << "\n";
        }

        StatusFile::addInfo("  Dropped " + to_string(lstats.MtargetOnly) + " variants not found in reference.<br>", false);
        yamlinfo << "    Dropped target only variants: " << lstats.MtargetOnly << "\n";
        if (lstats.MtargetOnly > M) {
            StatusFile::addWarning("More than 50% of the target variants not found in reference! Wrong genome built?");
        } else if (lstats.MtargetOnly > (M+lstats.MtargetOnly) / 10) {
            stringstream s;
            s << "More than 10% of the target variants not found in reference.";
            if (!allowRefAltSwap && !allowStrandFlip)
                s << " Allow ref/alt swaps and/or strand flips?";
            StatusFile::addWarning(s.str());
        }
        if (lstats.MmultiAllelicTgt) {
            StatusFile::addInfo("  Dropped " + to_string(lstats.MmultiAllelicTgt) + " multi-allelic variants in target.<br>", false);
            yamlinfo << "    Dropped multi-allelic target variants: " << lstats.MmultiAllelicTgt << "\n";
        }
        if (lstats.MmultiAllelicRefTgt) {
            if (excludeMultiAllRef) {
                StatusFile::addInfo("  Dropped " + to_string(lstats.MmultiAllelicRefTgt) + " variants bi-allelic in target but multi-allelic in reference.<br>", false);
                yamlinfo << "    Dropped multi-allelic reference variants: " << lstats.MmultiAllelicRefTgt << "\n";
            } else {
//                if (!loadQuickRef){
//                    StatusFile::addInfo("  Split " + to_string(lstats.MmultiAllelicRefTgt) + " variants bi-allelic in target but multi-allelic in reference into " + to_string(lstats.MmultiSplittedTgt)
//                            + " bi-allelic reference variants,<br>\n"
//                            + "  "  + to_string(lstats.GmultiFilledTgt) + " haplotypes filled with reference alleles.<br>", false);
//                    yamlinfo << "    Multi-allelic splits in common variants:\n";
//                    yamlinfo << "      Common multi-allelic reference variants before split: " << lstats.MmultiAllelicRefTgt << "\n";
//                    yamlinfo << "      Resulting bi-allelic reference variants: " << lstats.MmultiSplittedTgt << "\n";
//                    yamlinfo << "      Haplotypes filled with reference alleles: " << lstats.GmultiFilledTgt << "\n";
//                } else {
                StatusFile::addInfo("  " + to_string(lstats.MmultiAllelicRefTgt) + " variants bi-allelic in target were multi-allelic in reference.<br>", false);
                yamlinfo << "    Common reference variants from multi-allelic splits: " << lstats.MmultiAllelicRefTgt << "\n";
//                }
            }
        }
        if (lstats.MmonomorphicRefTgt) {
            StatusFile::addInfo("  Dropped " + to_string(lstats.MmonomorphicRefTgt) + " variants bi-allelic in target but monomorphic in reference.<br>", false);
            yamlinfo << "    Dropped target variants monomorphic in reference: " << lstats.MmonomorphicRefTgt << "\n";
        }
        if (lstats.MrefAltError) {
            StatusFile::addInfo("  Dropped " + to_string(lstats.MrefAltError) + " variants with allele mismatches.<br>", false);
            yamlinfo << "    Dropped allele mismatched variants: " << lstats.MrefAltError << "\n";
        }
        if (useExclude) {
            StatusFile::addInfo("  Dropped " + to_string(lstats.Mexclude) + " variants from target based on --vcfExclude.<br>", false);
            yamlinfo << "    User excluded variants: " << lstats.Mexclude << "\n";
        }

        if (lstats.Munphased) {
            StatusFile::addInfo("  Unphased variants in phasing output: " + to_string(lstats.Munphased) + "<br>", false);
            yamlinfo << "    Unphased variants after phasing: " << lstats.Munphased << "\n";
        }

        StatusFile::addInfo("</p><p class='pinfo'>", false);
//        size_t mrefonly = loadQuickRef ? (Mref-M) : lstats.MrefOnly;
        size_t mrefonly = Mref-M;
        StatusFile::addInfo("  " + to_string(mrefonly) + " variants in reference but not in target.<br>", false);
        yamlinfo << "    Reference-only variants: " << mrefonly << "\n";
        if (doImputation) {
            if (excludeMultiAllRef) {
                StatusFile::addInfo("  Excluding " + to_string(lstats.MmultiAllelicRefOnly) + " reference-only multi-allelic variants from imputation.<br>", false);
                yamlinfo << "    Excluded reference-only multi-allelic variants: " << lstats.MmultiAllelicRefOnly << "\n";
            } else if (lstats.MmultiAllelicRefOnly) {
//                if (!loadQuickRef){
//                    StatusFile::addInfo("  Split " + to_string(lstats.MmultiAllelicRefOnly) + " reference-only multi-allelic variants into "
//                            + to_string(lstats.MmultiSplittedRefOnly) + " bi-allelic reference variants,<br>\n"
//                            + "    " + to_string(lstats.GmultiFilledRefOnly) + " haplotypes filled with reference alleles.<br>", false);
//                    yamlinfo << "    Multi-allelic splits in reference-only variants:\n";
//                    yamlinfo << "      Reference-only multi-allelic variants before split: " << lstats.MmultiAllelicRefOnly << "\n";
//                    yamlinfo << "      Resulting bi-allelic reference variants: " << lstats.MmultiSplittedRefOnly << "\n";
//                    yamlinfo << "      Haplotypes filled with reference alleles: " << lstats.GmultiFilledRefOnly << "\n";
//                } else {
                StatusFile::addInfo("  " + to_string(lstats.MmultiAllelicRefOnly) + " bi-allelic variants in reference only that originate from multi-allelic splits.<br>", false);
                yamlinfo << "    Reference-only bi-allelic variants from multi-allelic splits: " << lstats.MmultiAllelicRefOnly << "\n";
//                }
            } else {
                StatusFile::addInfo("  No multi-allelic reference-only variants.<br>", false);
                // we skip this output in YAML
            }
        }

        if (M <= 1ull) {
            StatusFile::addWarning("<b>Analysis skipped:</b> Target and ref have too few matching variants. (M = " + to_string(M) + ")");
            StatusFile::updateStatus(1, "Skipped");
            if (yaml)
                StatusFile::addInfoYAML("Chunk", yamlinfo.str());
            // TODO This is not good! Instead, we should continue with the next chunk somehow, but what to do with the (reference) data covered by this chunk?
            exit(EXIT_SUCCESS);
        }

        StatusFile::addInfo("</p><p class='pinfo'>", false);

        if (lstats.MwithMissingRef) {
            StatusFile::addInfo("  Missing genotypes in phasing reference: " + to_string(lstats.GmissingRef) + " in " + to_string(lstats.MwithMissingRef) + " variants.<br>", false);
            yamlinfo << "    Missing genotypes in phasing reference:\n";
            yamlinfo << "      Number:" << lstats.GmissingRef << "\n";
            yamlinfo << "      Fraction: " << ((lstats.GmissingRef / (double) M) / Nref) << "\n";
            yamlinfo << "      Affected variants: " << lstats.MwithMissingRef << "\n";
            yamlinfo << "      Affected variants fraction: " << (lstats.MwithMissingRef / (double) M) << "\n";
            StatusFile::addWarning("Reference for phasing contains missing genotypes (set randomly according to allele frequency).<br>\n"
                    "   Variants with missing data: " + to_string(lstats.MwithMissingRef)
                    + " Fraction: " + to_string(lstats.MwithMissingRef / (double) M) + "<br>\n"
                    "   Missing genotypes: " + to_string(lstats.GmissingRef)
                    + " Fraction: " + to_string((lstats.GmissingRef / (double) M) / Nref));
        }
        if (lstats.MwithUnphasedRef) {
            StatusFile::addInfo("  Unphased genotypes in phasing reference: " + to_string(lstats.GunphasedRef) + " in " + to_string(lstats.MwithUnphasedRef) + " variants.<br>", false);
            yamlinfo << "    Unphased genotypes in phasing reference:\n";
            yamlinfo << "      Number: " << lstats.GunphasedRef << "\n";
            yamlinfo << "      Fraction: " << ((lstats.GunphasedRef / (double) M) / Nref) << "\n";
            yamlinfo << "      Affected variants: " << lstats.MwithUnphasedRef << "\n";
            yamlinfo << "      Affected variants fraction: " << (lstats.MwithUnphasedRef / (double) M) << "\n";
            StatusFile::addWarning("Reference for phasing contains unphased genotypes (set to random phase).<br>\n"
                    "   Variants with unphased data: " + to_string(lstats.MwithUnphasedRef)
                    + " Fraction: " + to_string(lstats.MwithUnphasedRef / (double) M) + "<br>\n"
                    "   Unphased genotypes: " + to_string(lstats.GunphasedRef)
                    + " Fraction: " + to_string((lstats.GunphasedRef / (double) M) / Nref));
        }

        if (doImputation) {
            if (lstats.MmissingRefOnly) {
                StatusFile::addInfo("  Missing genotypes in imputation reference: " + to_string(lstats.GmissingRefOnly) + " in " + to_string(lstats.MmissingRefOnly) + " variants.<br>", false);
                yamlinfo << "    Missing genotypes in imputation reference:\n";
                yamlinfo << "      Number: " << lstats.GmissingRefOnly << "\n";
                yamlinfo << "      Fraction: " << ((lstats.GmissingRefOnly / (double) Mref) / Nref) << "\n";
                yamlinfo << "      Affected variants: " << lstats.MmissingRefOnly << "\n";
                yamlinfo << "      Affected variants fraction: " << (lstats.MmissingRefOnly / (double) Mref) << "\n";
                StatusFile::addWarning("Reference for imputation contains missing genotypes (set randomly according to allele frequency).<br>\n"
                        "   Variants with missing data: " + to_string(lstats.MmissingRefOnly)
                        + " Fraction: " + to_string(lstats.MmissingRefOnly / (double) Mref) + "<br>\n"
                        "   Missing genotypes: " + to_string(lstats.GmissingRefOnly)
                        + " Fraction: " + to_string((lstats.GmissingRefOnly / (double) Mref) / Nref));
            }
            if (lstats.MunphasedRefOnly) {
                StatusFile::addInfo("  Unphased genotypes in imputation reference: " + to_string(lstats.GunphasedRefOnly) + " in " + to_string(lstats.MunphasedRefOnly) + " variants.<br>", false);
                yamlinfo << "    Unphased genotypes in imputation reference:\n";
                yamlinfo << "      Number: " << lstats.GunphasedRefOnly << "\n";
                yamlinfo << "      Fraction: " << ((lstats.GunphasedRefOnly / (double) Mref) / Nref) << "\n";
                yamlinfo << "      Affected variants: " << lstats.MunphasedRefOnly << "\n";
                yamlinfo << "      Affected variants fraction: " << (lstats.MunphasedRefOnly / (double) Mref) << "\n";
                StatusFile::addWarning("Reference for imputation contains unphased genotypes (set to random phase).<br>\n"
                        "   Variants with unphased data: " + to_string(lstats.MunphasedRefOnly)
                        + " Fraction: " + to_string(lstats.MunphasedRefOnly / (double) Mref) + "<br>\n"
                        "   Unphased genotypes: " + to_string(lstats.GunphasedRefOnly)
                        + " Fraction: " + to_string((lstats.GunphasedRefOnly / (double) Mref) / Nref));
            }
        }

        double missingrate = (lstats.GmissingTarget / (double) M) / Ntarget;
        StatusFile::addInfo("Av. missing rate in target genotypes: " + to_string(missingrate) + "</p>", false);
        yamlinfo << "    Missing rate in target genotypes: " << missingrate << "\n";
        if (missingrate > 0.1)
            StatusFile::addWarning("Missing rate in target genotypes is >10%.");

        if (Mref < 3000)
            StatusFile::addWarning("Less than 3000 reference variants in chunk.");

        // write YAML
        if (yaml)
            StatusFile::addInfoYAML("Chunk", yamlinfo.str());

        // collect local statistics for global statistics
        globalstats += lstats;

    } else { // create QRef


//        StatusFile::updateStatus(0, "write Qref meta");

        qrefvarsofs.close();
        qRefWriteMetaAndConcat(); // dump the reference meta data to qref file and concat with variant data to create final qref

        // statistics for qref creation
        cout << "Variants dumped: " << Mref << endl;

        if (lstats.MmultiAllelicRefOnly) {
            StatusFile::addWarning("Split " + to_string(lstats.MmultiAllelicRefOnly)
                    + " reference-only multi-allelic SNPs into " + to_string(lstats.MmultiSplittedRefOnly)
                    + " bi-allelic reference SNPs,<br>\n"
                    "  " + to_string(lstats.GmultiFilledRefOnly) + " haplotypes filled with reference alleles.");
        }
        if (lstats.MmissingRefOnly) {
            StatusFile::addWarning("Reference for imputation contains missing genotypes (set randomly according to allele frequency).<br>\n"
                    "   Sites with missing data: " + to_string(lstats.MmissingRefOnly)
                    + " Fraction: " + to_string(lstats.MmissingRefOnly / (double) Mref) + "<br>\n"
                    "   Missing genotypes: " + to_string(lstats.GmissingRefOnly)
                    + " Fraction: " + to_string(lstats.GmissingRefOnly / (double) Mref / Nref));
        }
        if (lstats.MunphasedRefOnly) {
            StatusFile::addWarning("Reference for imputation contains unphased genotypes (set to random phase).<br>\n"
                    "   Sites with unphased data: " + to_string(lstats.MunphasedRefOnly)
                    + " Fraction: " + to_string(lstats.MunphasedRefOnly / (double) Mref) + "<br>\n"
                    "   Unphased genotypes: " + to_string(lstats.GunphasedRefOnly)
                    + " Fraction: " + to_string(lstats.GunphasedRefOnly / (double) Mref / Nref));
        }
    }

    StatusFile::nextStep();
}

// keep SNP for imputation
inline void VCFData::processReferenceOnlySNP(bcf1_t *ref, void **ref_gt, int *mref_gt, int &numMissing, int &numUnphased, bool multiallelic) {
    float af;
//    int nref_gt = bcf_get_genotypes(ref_hdr, ref, ref_gt, mref_gt);
    processReferenceSNP(Nref, ref, ref_gt, mref_gt, numMissing, numUnphased, af, false);
//    cout << "SNP: " << Mref << " Nref: " << Nref << " nref_gt: " << nref_gt << " nMiss: " << numMissing << " nUnph: " << numUnphased << endl;
    alleleFreqsFullRefRegion.push_back(af);
    // store position information of variant
    positionsFullRefRegion.push_back(ref->pos);
    // store allele information of variant
    allelesFullRefRegion.emplace_back(ref->d.allele[0]);
    allelesFullRefRegion.emplace_back(ref->d.allele[1]);
    variantIDsFullRefRegion.emplace_back(ref->d.id);
    multiAllFlagsFullRefRegion.push_back_withPreInit(multiallelic);
    if (multiallelic)
        MrefMultiAllreg++;

    currMref++;

    // not required anymore as this function is called only when creating a Qref
//    isImputed.push_back(true);
//    nextPidx.push_back(bcf_pout.size()); // will point to the next common index or the end, if there will be no more common sites
//    ptgtIdx.push_back(currM);
//    if (inOverlap) {
//        isImputedOverlap.push_back(true);
//        nextPidxOverlap.push_back(bcf_poutOverlap.size()); // will point to the next common index or the end, if there will be no more common sites
//        ptgtIdxOverlap.push_back(currMOverlap);
//        currMrefOverlap++;
//    }
}

inline void VCFData::processReferenceSNP(int nsmpl, bcf1_t *ref, void **ref_gt, int *mref_gt,
		int &numMissing, int &numUnphased, float &af, bool addForPhasing) {
    numMissing = numUnphased = 0;

    // reserve space for haplotype data
    BooleanVector::data_type* hapdata = (BooleanVector::data_type*) MyMalloc::calloc(hapcapacity, 1, "hapdata_processRefSNP"); // need this pre-initialized with zeros!
    if (!hapdata) {
        cout << endl;
        StatusFile::addError(string("Not enough memory for reference. Alloc size: ").append(to_string(hapcapacity)));
        exit(EXIT_FAILURE);
    }
    referenceFullT.emplace_back(hapdata, hapcapacity, Nrefhaps);

    // set size of current variant in transposed reference to the number of haps to be added
    if (addForPhasing) {
        referenceT[currM].setSize(2*nsmpl);
    }

    // fetch haplotypes
    int ngt = bcf_get_genotypes(ref_hdr, ref, ref_gt, mref_gt);
    const int32_t *gt = reinterpret_cast<const int32_t*>(*ref_gt);
    if (ngt != 2 * nsmpl && ngt != nsmpl) {
        cout << endl;
        stringstream ss;
        ss << "Unallowed encoding. Reference encoding ploidy != 2 (ngt != 2*nsmpl) and != 1 (ngt != nsmpl): ngt="
           << ngt << ", nsmpl=" << nsmpl;
        StatusFile::addError(ss.str());
        exit(EXIT_FAILURE);
    }
    int ploidy = ngt == nsmpl ? 1 : 2; // others are not possible (see trap above)
    const int32_t *ptr = gt;
    int ac1 = 0, an = 0;
    bool refafvalid = false;
    float refaf = 0.0;
    for (int i = 0; i < nsmpl; i++, ptr += ploidy) {
        Haplotype haps[2] = {Haplotype::Ref, Haplotype::Ref};
        bool unphased = false;
        if (*ptr == bcf_int32_vector_end) {
            cout << endl;
            StatusFile::addError("Zero ploidy in reference?!");
            exit(EXIT_FAILURE);
        }
        if (ploidy == 2) { // normal case
            if (bcf_gt_is_missing(*ptr) || bcf_gt_is_missing(*(ptr+1))) { // missing allele
                // get the allele frequency directly from the reference (if not already done before)
                if (!refafvalid) {
                    // whatever happens here, we don't try a second time for this variant
                    refafvalid = true;
                    // you need to manually alloc the memory for the return value on the heap,
                    // as HTSlib may re-allocate the memory (e.g. if the AF tag has more than one entry (multi-allelics))
                    int refaf_size = sizeof(float);
                    float *refaf_ptr = (float*) malloc(refaf_size);
                    int numaf = bcf_get_info_float(ref_hdr, ref, "AF", (void*)&refaf_ptr, &refaf_size);
                    if (numaf <= 0) // no AF tag or error fetching the AF info
                        refaf = 0.5; // set to "fifty-fifty"
                    else
                        refaf = *refaf_ptr; // if there's more than one AF entry, take the first. This is ok, as it's only to fill the missing allele.
                    free(refaf_ptr);
                }
                numMissing++; // only counted per genotype
            }

            // first allele
            if (bcf_gt_is_missing(*ptr)) {
                // set missing haplotype randomly according to allele frequency
                if (deterministic)
                    haps[0] = refaf <= 0.5 ? Haplotype::Ref : Haplotype::Alt;
                else
                    haps[0] = randDist(randGen) > refaf ? Haplotype::Ref : Haplotype::Alt;
            } else { // not missing
                haps[0] = (bcf_gt_allele(*ptr) >= 1) ? Haplotype::Alt : Haplotype::Ref; // encode REF allele -> 0, ALT allele(s) -> 1
            }
            // count allele
            if (haps[0] == Haplotype::Alt)
                ac1++;
            an++;

            // second allele
            if (bcf_gt_is_missing(*(ptr+1))) {
                if (deterministic)
                    haps[1] = refaf <= 0.5 ? Haplotype::Ref : Haplotype::Alt;
                else
                    haps[1] = randDist(randGen) > refaf ? Haplotype::Ref : Haplotype::Alt;
                // count allele
                if (haps[1] == Haplotype::Alt)
					ac1++;
                an++;
            } else if (*(ptr+1) == bcf_int32_vector_end) {
                // haploid or missing
//                if (allowHaploid || bcf_gt_is_missing(*ptr)) {
                    // simply encode as homozygous diploid
                    haps[1] = haps[0];
                    if (!bcf_gt_is_missing(*ptr)) { // not missing -> haploid sample
                        if (haploidsRef_initialized[i]) {
                            if (!haploidsRef[i]) { // already initialized as diploid
                                cout << endl;
                                StatusFile::addError("Reference sample " + to_string(i) + " was diploid and is haploid now! Incorrect PAR regions?");
                                exit(EXIT_FAILURE);
                            }
                        } else { // not initialized yet
                            haploidsRef_initialized[i] = true;
                            haploidsRef[i] = true;
                            nHaploidsRef++;
                        }
                    }
//                } else {
//                    cout << endl;
//                    StatusFile::addError("Reference contains haploid sample.");
//                    exit(EXIT_FAILURE);
//                }
            } else { // not missing and not haploid -> diploid
                haps[1] = (bcf_gt_allele(*(ptr+1)) >= 1) ? Haplotype::Alt : Haplotype::Ref; // encode REF allele -> 0, ALT allele(s) -> 1
                if (!bcf_gt_is_phased(*(ptr+1)))
                    unphased = true;
//                if (allowHaploid) {
                    if (haploidsRef_initialized[i]) {
                        if (haploidsRef[i]) {
                            cout << endl;
                            StatusFile::addError("Reference sample " + to_string(i) + " was haploid and is diploid now! Incorrect PAR regions?");
                            exit(EXIT_FAILURE);
                        }
                    } else {
                        // already pre-initialized as diploid... haploidsRef[i] = false; // diploid
                        haploidsRef_initialized[i] = true;
                    }
//                }
				// count allele
				if (haps[1] == Haplotype::Alt)
					ac1++;
				an++;
            }

        } else { // ploidy == 1
            if (bcf_gt_is_missing(*ptr)) { // missing allele
                // get the allele frequency directly from the reference (if not already done before)
                if (!refafvalid) {
                    // whatever happens here, we don't try a second time for this variant
                    refafvalid = true;
                    // you need to manually alloc the memory for the return value on the heap,
                    // as HTSlib may re-allocate the memory (e.g. if the AF tag has more than one entry (multi-allelics))
                    int refaf_size = sizeof(float);
                    float *refaf_ptr = (float*) malloc(refaf_size);
                    int numaf = bcf_get_info_float(ref_hdr, ref, "AF", (void*)&refaf_ptr, &refaf_size);
                    if (numaf <= 0) // no AF tag or error fetching the AF info
                        refaf = 0.5; // set to "fifty-fifty"
                    else
                        refaf = *refaf_ptr; // if there's more than one AF entry, take the first. This is ok, as it's only to fill the missing allele.
                    free(refaf_ptr);
                }
                // set missing haplotype randomly according to allele frequency (and encode as homozygous diploid)
                if (deterministic)
                    haps[0] = haps[1] = refaf <= 0.5 ? Haplotype::Ref : Haplotype::Alt;
                else
                    haps[0] = haps[1] = randDist(randGen) > refaf ? Haplotype::Ref : Haplotype::Alt;
                numMissing++;
            } else {
                haps[0] = haps[1] = (bcf_gt_allele(*ptr) >= 1) ? Haplotype::Alt : Haplotype::Ref; // encode as diploid homozygous
                if (haploidsRef_initialized[i]) {
                    if (!haploidsRef[i]) {
                        cout << endl;
                        StatusFile::addError("Reference sample " + to_string(i) + " was diploid and is haploid now! Incorrect PAR regions?");
                        exit(EXIT_FAILURE);
                    }
                } else {
                    haploidsRef[i] = true; // haploid
                    haploidsRef_initialized[i] = true;
                    nHaploidsRef++;
                }
            }
            // count allele
            if (haps[0] == Haplotype::Alt)
                ac1++;
            an++;
        }

        if (unphased) {
            if (haps[0] != haps[1] && !deterministic && randDist(randGen) > 0.5)
                swap(haps[0], haps[1]); // randomize phasing
            numUnphased++;
        }
        // add to vectors
        if (addForPhasing) {
            reference[2*i].push_back_withPreInit(toBool(haps[0]));
            reference[2*i+1].push_back_withPreInit(toBool(haps[1]));
            referenceT[currM].setPairWithPreInit(2*i, toBool(haps[0]), toBool(haps[1]));
        }
        if (doImputation || createQRef) {
            referenceFullT.back().setPairWithPreInit(2*i, toBool(haps[0]), toBool(haps[1]));
        }
    }

    // compute allele frequency
    af = 0;
    if (an)
        af = ((float) ac1) / an;

    // add to Qref
    if (createQRef) {
        qRefAppendVariant(referenceFullT.back());
        // we don't need the data in memory anymore if we just create the Qref
        referenceFullT.clear();
        MyMalloc::free(hapdata);
    }
}

inline void VCFData::processTargetSNP(int nsmpl, int ngt, const int32_t *gt, bool refAltSwap, int &numMissing, size_t tgtvariantpos) {
    numMissing = 0;
    if (ngt != 2 * nsmpl && ngt != nsmpl) {
        cout << endl;
        stringstream ss;
        ss << "Unallowed encoding. Target encoding ploidy != 2 (ngt != 2*nsmpl) and != 1 (ngt != nsmpl): ngt="
           << ngt << ", nsmpl=" << nsmpl;
        StatusFile::addError(ss.str());
        exit(EXIT_FAILURE);
    }

    int ploidy = ngt == nsmpl ? 1 : 2; // others are not possible (see trap above)
    const int32_t *ptr = gt;
    for (int i = 0; i < nsmpl; i++, ptr += ploidy) {
        Genotype g = Genotype::HomRef;
        bool phase = false;
        if (*ptr == bcf_int32_vector_end) {
            cout << endl;
            StatusFile::addError("Zero ploidy in target?!");
            exit(EXIT_FAILURE);
        }
        if (ploidy == 2) { // normal case
            if (bcf_gt_is_missing(*ptr) || bcf_gt_is_missing(*(ptr+1))) { // missing allele
                g = Genotype::Miss;
                if (noImpMissing) { // will only be required if we are not going to impute missings during phasing
                    tgtmiss[i].push_back(targets[i].size()); // store position of missing genotype
                }
                numMissing++;
            } else if (*(ptr+1) == bcf_int32_vector_end) { // haploid sample!
//                if (allowHaploid) { // on chrX/Y: encode as homozygous diploid and mark as male
                    g = bcf_gt_allele(*ptr) >= 1 ? Genotype::HomAlt : Genotype::HomRef;
                    // mark as haploid
                    if (haploidsTgt_initialized[i]) {
                        if (!haploidsTgt[i]) {
                            cout << endl;
                            stringstream ss;
                            ss << "Target sample " << i << " (" << targetIDs[i] << " at " << tgtvariantpos+1 << ") was diploid and is haploid now! Incorrect PAR regions?";
                            StatusFile::addError(ss.str());
                            exit(EXIT_FAILURE);
                        }
                    } else {
                        haploidsTgt[i] = true; // haploid
                        haploidsTgt_initialized[i] = true;
                        nHaploidsTgt++;
                        if (iters > 1) {
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                            if (i >= DEBUG_TARGET_START && i <= DEBUG_TARGET_STOP)
                                haploidsRef[Nref+i-DEBUG_TARGET_START] = true; // this target will be ref in a phasing iteration > 1
#else
                            haploidsRef[Nref+i] = true; // this target will be ref in a phasing iteration > 1
#endif
                        }
                    }
//                } else {
//                    cout << endl;
//                    StatusFile::addError("Target contains haploid sample.");
//                    exit(EXIT_FAILURE);
//                }
            } else { // diploid
                int idx =  bcf_gt_allele(*ptr) + bcf_gt_allele(*(ptr+1));
                switch (idx) {
                case 0:
                    g = Genotype::HomRef;
                    break;
                case 1:
                    g = Genotype::Het;
                    phase = bcf_gt_allele(*ptr) != 0;
                    break;
                default: // no other possibility (multi-allelic cases have been filtered before)
                    g = Genotype::HomAlt;
                }
//                if (allowHaploid) {
                    if (haploidsTgt_initialized[i]) {
                        if (haploidsTgt[i]) {
                            cout << endl;
                            stringstream ss;
                            ss << "Target sample " << i << " (" << targetIDs[i] << " at " << tgtvariantpos+1 << ") was diploid and is haploid now! Incorrect PAR regions?";
                            StatusFile::addError(ss.str());
                            exit(EXIT_FAILURE);
                        }
                    } else {
                        // already pre-initialized as diploid... haploidsTgt[i] = false; // diploid
                        haploidsTgt_initialized[i] = true;
                    }
//                }
            }
        } else { // ploidy == 1
            if (bcf_gt_is_missing(*ptr)) { // missing allele
                g = Genotype::Miss;
                if (noImpMissing) { // will only be required if we are not going to impute missings during phasing
                    tgtmiss[i].push_back(targets[i].size()); // store position of missing genotype
                }
                numMissing++;
            } else { // encode as homozygous diploid and mark as haploid
                g = bcf_gt_allele(*ptr) >= 1 ? Genotype::HomAlt : Genotype::HomRef;
                if (haploidsTgt_initialized[i]) {
                    if (!haploidsTgt[i]) {
                        cout << endl;
                        stringstream ss;
                        ss << "Target sample " << i << " (" << targetIDs[i] << " at " << tgtvariantpos+1 << ") was diploid and is haploid now! Incorrect PAR regions?";
                        StatusFile::addError(ss.str());
                        exit(EXIT_FAILURE);
                    }
                } else {
                    haploidsTgt[i] = true; // haploid
                    haploidsTgt_initialized[i] = true;
                    nHaploidsTgt++;
                }
            }
        }

        if (refAltSwap) { // the target's alleles are swapped relative to the reference ones, and we stick to the reference!
            switch (g) {
            case Genotype::HomRef:
                g = Genotype::HomAlt;
                break;
            case Genotype::HomAlt:
                g = Genotype::HomRef;
                break;
            default:
                break;
            }
        }

        targets[i].push_back(g);
        if (skipPhasing) { // store phase
            tgtinphase[i].push_back(phase);
        }
    }
}

inline void VCFData::qRefWriteMetaAndConcat() {
    // dump the reference meta data to qref meta file, then concat qref variant file

    ofstream qout(qreffilename, ios_base::binary);
    // header
    qout << "QREF"; // "magic"
    char tmp[4];
    tmp[0] = (char) QREFV_MAJ; // major version
    tmp[1] = (char) QREFV_MIN; // minor version
    tmp[2] = 0;
    tmp[3] = (char) chrom; // chromosome
    qout.write(tmp,4);
    qout.write((char*)&Nref, sizeof(size_t));
    qout.write((char*)&nHaploidsRef, sizeof(size_t));
    qout.write((char*)&Mref, sizeof(size_t));
    qout.write((char*)&MrefMultiAllreg, sizeof(size_t));
    qout << flush;

    // if this is chrX, dump the information on haploid reference samples (for Y we know that all samples are haploid)
    // Note, this is kept to be backward compatible.
    // In general, if there is at least one haploid sample, we dump the flags, with the only exception if ALL samples are haploid (e.g. for chrY).
    if ((nHaploidsRef > 0 && nHaploidsRef != Nref) || chrom == CHRX) {
        // we spend one byte per flag... TODO could be compacted!
        for (bool flag : haploidsRef) {
            qout << (char)(flag ? -1 : 0);
        }
        if (qout.fail()) {
            StatusFile::addError("Write error while writing haploid flags.");
            exit(EXIT_FAILURE);
        }
    }

    // variant positions
    qout.write((char*)positionsFullRefRegion.data(), Mref*sizeof(int64_t));
    if (qout.fail()) {
        StatusFile::addError("Write error while writing variant positions.");
        exit(EXIT_FAILURE);
    }

    // allele frequencies
    qout.write((char*)alleleFreqsFullRefRegion.data(), Mref*sizeof(float));
    if (qout.fail()) {
        StatusFile::addError("Write error while writing allele frequencies.");
        exit(EXIT_FAILURE);
    }

    // allele strings
    for (const string &all: allelesFullRefRegion) {
        qout.write(all.c_str(), all.length()+1); // +1: don't forget the zero termination character!
        if (qout.fail()) {
            StatusFile::addError("Write error while writing allele strings.");
            exit(EXIT_FAILURE);
        }
    }

    // allele IDs
    for (const string &all: variantIDsFullRefRegion) {
        qout.write(all.c_str(), all.length()+1); // +1: don't forget the zero termination character!
        if (qout.fail()) {
            StatusFile::addError("Write error while writing allele IDs.");
            exit(EXIT_FAILURE);
        }
    }

    // flags for multi-allelic splits
    size_t macapacity = roundToMultiple<size_t>(Mref, UNITWORDS*sizeof(BooleanVector::data_type)*8)/8; // space for Mref SNPs
    qout.write((char*)multiAllFlagsFullRefRegion.getData(), macapacity);
    if (qout.fail()) {
        StatusFile::addError("Write error while writing multi-allelic flags.");
        exit(EXIT_FAILURE);
    }

    qout.close();

    // concat with variant data
    vector<string> concatfiles(2);
    concatfiles[0] = qreffilename;
    concatfiles[1] = qreffilename+".vars";
    concatFiles(concatfiles);
}

inline void VCFData::qRefAppendVariant(const BooleanVector& var) {

    // haplotype data (only the final Mref ones! size could be larger due to reservations for multi-allelic splits!)
    vector<char> rlenc; // destination vector for runlength encodings
    bool encode_success = runlengthEncode(var.getData(), hapcapacity/sizeof(BooleanVector::data_type), rlenc);
    // if encoding was not successful, the length is encoded with zero to indicate that the original sequence will be written
    size_t rlencsize = encode_success ? rlenc.size() : 0ull; // if encoding was not successful, the length is encoded with zero to indicate that the original sequence will be written
    qrefvarsofs.write((const char*)&rlencsize, sizeof(size_t)); // write the size of the vector
    if (encode_success)
        qrefvarsofs.write(rlenc.data(), rlenc.size()); // write the encoded vector itself after successful encoding
    else
        qrefvarsofs.write((const char*)(var.getData()), hapcapacity); // write the original sequence if runlength encoding would lead to a longer sequence
    if (qrefvarsofs.fail()) {
        cout << endl;
        StatusFile::addError("Write error while writing haplotypes.");
        exit(EXIT_FAILURE);
    }

}

inline void VCFData::qRefOpenReadMeta() {
//    StatusFile::updateStatus(0, "Loading Qref metadata");

    qin.open(refFileName, ios_base::binary);
    if (qin.fail()) {
        StatusFile::addError("Could not open Qref file.");
        exit(EXIT_FAILURE);
    }
    // read header
    char head[8];
    qin.read(head, 8);
    // check header
    if (head[0] != 'Q' || head[1] != 'R' || head[2] != 'E' || head[3] != 'F' || head[6] != 0) {
        StatusFile::addError("Not a valid Qref file.");
        exit(EXIT_FAILURE);
    }
    if (head[4] != (char) QREFV_MAJ || head[5] > (char) QREFV_MIN) { // major version must be equal, minor must not be newer!
        StatusFile::addError(string("Wrong Qref version. Expected v") + to_string(QREFV_MAJ) + "." + to_string(QREFV_MIN) + ".");
        exit(EXIT_FAILURE);
    }

    // now we believe this is a valid qref file

    // chromosome, number of samples, number of SNPs
    chrom = (int) head[7]; // TODO Note: chromlit is taken from the filename yet... change this here?
    qin.read((char*)&Nref, sizeof(size_t));
    qin.read((char*)&nHaploidsRef, sizeof(size_t));
    qin.read((char*)&Mrefglob, sizeof(size_t));
    qin.read((char*)&MrefMultiAllglob, sizeof(size_t));
    Nrefhaps = 2*Nref;
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
    Nrefhapsmax = iters > 1 ? Nrefhaps + 2*(DEBUG_TARGET_STOP-DEBUG_TARGET_START+1) : Nrefhaps;
#else
    Nrefhapsmax = iters > 1 ? Nrefhaps + 2*Ntarget : Nrefhaps;
#endif
    Mrefreg = Mrefglob; // will be updated below if we have selected a region

    if (Mrefglob == 0) // stop, if there's no data
        return;

    // if this is chrX, load the information on haploid reference samples (for Y and others we know that all samples are haploid/diploid)
    // Note, we only load this information per default from chrX to be backward compatible. Otherwise, we load this if at least one sample
    // is haploid with the exception if we know that ALL samples are haploid.
    haploidsRef.clear();
    haploidsRef.resize(Nrefhapsmax/2, nHaploidsRef == Nref); // initialized with all diploid if not all samples are haploid
    haploidsRef_initialized.clear();
    haploidsRef_initialized.resize(Nref, true); // all initialized, since we are loading them from Qref
    haploidsRefMap.clear();
    haploidsRefMap.reserve(Nrefhapsmax);
    if ((nHaploidsRef > 0 && nHaploidsRef != Nref) || chrom == CHRX) {
        for (size_t i = 0; i < Nref; i++) {
            char flag;
            qin >> flag;
            if (qin.fail() || qin.eof()) {
                StatusFile::addError("Failed reading haploid flags from reference.");
                exit(EXIT_FAILURE);
            }
            haploidsRef[i] = flag != 0;
            haploidsRefMap.push_back(2*i);
            if (flag == 0) // add another index only if not haploid
                haploidsRefMap.push_back(2*i+1);
        }
        // NOTE: the missing flags (Nrefhaps to Nrefhapsmax) are created during reading of target!
    } else if (nHaploidsRef == Nref) {
        // all haploid: map to every second index
        for (size_t i = 0; i < Nrefhapsmax/2; i++)
            haploidsRefMap.push_back(2*i);
    } else { // all diploid chromosome: identity
        for (size_t i = 0; i < Nrefhapsmax; i++)
            haploidsRefMap.push_back(i);
    }

    // load variant positions -> determine indices to load to fulfill region requests
    vector<int64_t> positionsFullRefpre(Mrefglob);
    size_t readbytes = 0;
    while (readbytes < Mrefglob*sizeof(int64_t)) {
        qin.read(((char*)positionsFullRefpre.data())+readbytes, Mrefglob*sizeof(int64_t) - readbytes);
        if (qin.fail() || qin.eof()) {
            StatusFile::addError("Failed reading variant positions from reference.");
            exit(EXIT_FAILURE);
        }
        readbytes += qin.gcount();
    }

    // find indices
    qrefregfirstidx = 0;
    while (positionsFullRefpre[qrefregfirstidx] < startFlankRegionBp)
        qrefregfirstidx++;
    qrefreglastidx = Mrefglob-1;
    while (positionsFullRefpre[qrefreglastidx] > endFlankRegionBp)
        qrefreglastidx--;
    qrefreglastidx++; // is exclusive now
    Mrefreg = qrefreglastidx - qrefregfirstidx;
    if (Mrefreg == 0) // stop, if there's no data
        return;
    if (Mrefreg == Mrefglob) // all variants are in the region
        positionsFullRefRegion = move(positionsFullRefpre);
    else // copy only subset
        positionsFullRefRegion.assign(positionsFullRefpre.begin()+qrefregfirstidx, positionsFullRefpre.begin()+qrefreglastidx);

    // allele frequencies
    alleleFreqsFullRefRegion.clear();
    alleleFreqsFullRefRegion.resize(Mrefreg);
    qin.seekg(qrefregfirstidx*sizeof(float), ios_base::cur); // skip the first fields until region
    readbytes = 0;
    while (readbytes < Mrefreg*sizeof(float)) {
        qin.read(((char*)alleleFreqsFullRefRegion.data())+readbytes, Mrefreg*sizeof(float) - readbytes);
        if (qin.fail() || qin.eof()) {
            StatusFile::addError("Failed reading allele frequencies from reference.");
            exit(EXIT_FAILURE);
        }
        readbytes += qin.gcount();
    }
    qin.seekg((Mrefglob-qrefreglastidx)*sizeof(float), ios_base::cur); // skip the fields after region

    // allele strings
    allelesFullRefRegion.clear();
    allelesFullRefRegion.reserve(2*Mrefreg);
    // ignore the first strings until region
    for (size_t i = 0; i < 2*qrefregfirstidx; i++)
        qin.ignore(numeric_limits<streamsize>::max(), 0); // delimiting character is the zero termination char
    char tmpbuf[TMPBUFSIZE]; // temporary buffer
    tmpbuf[TMPBUFSIZE-1] = 0; // ensure last char is the zero termination char
    bool appendflag = false;
    while (qin && allelesFullRefRegion.size() < 2*Mrefreg) {
        qin.getline(tmpbuf, TMPBUFSIZE, 0); // zero byte is the delimiting character
        if (qin.fail() && qin.gcount() != TMPBUFSIZE-1) { // real failure
            StatusFile::addError("Failed reading allele strings from reference.");
            exit(EXIT_FAILURE);
        }
        // else if failbit is set, the buffer was not large enough, so exactly TMPBUFSIZE-1 characters were read and can be normally processed (buf[TMPBUFSIZE-1]==0)
        if (appendflag)
            allelesFullRefRegion.back().append(tmpbuf);
        else
            allelesFullRefRegion.emplace_back(tmpbuf);
        appendflag = false;
        if (qin.fail() && qin.gcount() == TMPBUFSIZE-1) { // buffer was not large enough -> append after next read and clear failbit
            appendflag = true;
            qin.clear();
        }
    }
    if (allelesFullRefRegion.size() != 2*Mrefreg) {
        StatusFile::addError("Failed reading allele strings from reference (size != 2*Mrefreg).");
        exit(EXIT_FAILURE);
    }
    // ignore strings after region
    for (size_t i = 2*qrefreglastidx; i < 2*Mrefglob; i++)
        qin.ignore(numeric_limits<streamsize>::max(), 0); // delimiting character is the zero termination char

    // allele IDs
    variantIDsFullRefRegion.clear();
    variantIDsFullRefRegion.reserve(Mrefreg);
    // ignore the first strings until region
    for (size_t i = 0; i < qrefregfirstidx; i++)
        qin.ignore(numeric_limits<streamsize>::max(), 0); // delimiting character is the zero termination char
    appendflag = false;
    while (qin && variantIDsFullRefRegion.size() < Mrefreg) {
        qin.getline(tmpbuf, TMPBUFSIZE, 0); // zero byte is the delimiting character
        if (qin.fail() && qin.gcount() != TMPBUFSIZE-1) { // real failure
            StatusFile::addError("Failed reading allele IDs from reference.");
            exit(EXIT_FAILURE);
        }
        // else if failbit is set, the buffer was not large enough, so exactly TMPBUFSIZE-1 characters were read and can be normally processed (buf[TMPBUFSIZE-1]==0)
        if (appendflag)
            variantIDsFullRefRegion.back().append(tmpbuf);
        else
            variantIDsFullRefRegion.emplace_back(tmpbuf);
        appendflag = false;
        if (qin.fail() && qin.gcount() == TMPBUFSIZE-1) { // buffer was not large enough -> append after next read and clear failbit
            appendflag = true;
            qin.clear();
        }
    }
    if (variantIDsFullRefRegion.size() != Mrefreg) {
        StatusFile::addError("Failed reading allele IDs from reference (size != Mrefreg).");
        exit(EXIT_FAILURE);
    }

    // ignore strings after region
    for (size_t i = qrefreglastidx; i < Mrefglob; i++)
        qin.ignore(numeric_limits<streamsize>::max(), 0); // delimiting character is the zero termination char

    // flags for multi-allelic splits
    size_t macapacitypre = roundToMultiple<size_t>(Mrefglob, UNITWORDS*sizeof(BooleanVector::data_type)*8)/8; // space for Mrefglob SNPs
    BooleanVector::data_type* madatapre = (BooleanVector::data_type*) MyMalloc::malloc(macapacitypre, "madatapre"); // no pre-initialization required, contents will be overwritten

    readbytes = 0;
    while (readbytes < macapacitypre) {
        qin.read(((char*)madatapre)+readbytes, macapacitypre - readbytes);
        if (qin.fail() || qin.eof()) {
            StatusFile::addError("Failed reading multi-allelic flags from reference.");
            exit(EXIT_FAILURE);
        }
        readbytes += qin.gcount();
    }

    if (excludeMultiAllRef) { // only set the data if we need it
        if (Mrefglob == Mrefreg) { // no region selected -> take everything
            multiAllFlagsFullRefRegion.setData(madatapre, macapacitypre, Mrefreg);
            MrefMultiAllreg = MrefMultiAllglob;
        } else { // copy only selected region
            size_t macapacity = roundToMultiple<size_t>(Mrefreg, UNITWORDS*sizeof(BooleanVector::data_type)*8)/8; // space for Mrefreg SNPs
            BooleanVector::data_type* madata = (BooleanVector::data_type*) MyMalloc::malloc(macapacity, "madata_qrefOpen");
            BooleanVector maflagspre(madatapre, macapacitypre, Mrefglob);
            multiAllFlagsFullRefRegion.setDataAndInit(madata, macapacity, 0, false);
            MrefMultiAllreg = 0;
            for (size_t i = qrefregfirstidx; i < qrefreglastidx; i++) {
                multiAllFlagsFullRefRegion.push_back_withPreInit(maflagspre[i]);
                if (maflagspre[i])
                    MrefMultiAllreg++;
            }
        }
    } else { // multi-allelic reference variants are not excluded, so we set everything to false
        // ATTENTION! If a region was selected, the capacity does not necessarily conform with the size
        multiAllFlagsFullRefRegion.setDataAndInit(madatapre, macapacitypre, Mrefreg, false);
        MrefMultiAllreg = 0;
    }

    // now, the datastream points to the beginning of the samples' haplotypes
    hapcapacity = roundToMultiple<size_t>(Nrefhaps, UNITWORDS*sizeof(BooleanVector::data_type)*8)/8; // space for 2*Nref haps

    qrefcurridxglob = 0;
    // forward to first haps in selected region
    while(qrefcurridxglob < qrefregfirstidx) { // skip the first haps until beginning of region
        // read size of encoded data
        size_t encsize;
        qin.read((char*)&encsize, sizeof(size_t));
        // if the data was runlength encoded, encsize is > 0, else we need to skip the uncoded sequence
        if (encsize == 0)
            encsize = hapcapacity;
        // skip sequence
        qin.seekg(encsize, ios_base::cur);
        qrefcurridxglob++;
    }
    qrefcurridxreg = 0;
}


inline void VCFData::qRefLoadNextVariant() {
    // load haplotype data

    // reserve space for data
    BooleanVector::data_type* hapdata = (BooleanVector::data_type*) MyMalloc::calloc(hapcapacity, 1, "hapdata_qref"); // pre-initialization with zeroes
    if (!hapdata) {
        cout << endl;
        StatusFile::addError(string("Not enough memory for reference. Alloc size: ").append(to_string(hapcapacity)));
        exit(EXIT_FAILURE);
    }
    referenceFullT.emplace_back(hapdata, hapcapacity, Nrefhaps);

    // read size of encoded data
    size_t encsize;
    qin.read((char*)&encsize, sizeof(size_t));

    size_t readbytes = 0;
    if (encsize) { // encoding was successful and we read a runlength encoded sequence here
        // read encoded data
        vector<char> enc(encsize);
        while (readbytes < encsize) {
            qin.read(enc.data()+readbytes, encsize - readbytes);
            if (qin.fail() || qin.eof()) {
                cout << endl;
                StatusFile::addError("Failed reading haplotypes from reference.");
                exit(EXIT_FAILURE);
            }
            readbytes += qin.gcount();
        }
        // decode
        runlengthDecode(enc, hapdata, hapcapacity/sizeof(BooleanVector::data_type));
    } else { // the original sequence was stored
        while (readbytes < hapcapacity) {
            qin.read(((char*)hapdata)+readbytes, hapcapacity - readbytes);
            if (qin.fail() || qin.eof()) {
                cout << endl;
                StatusFile::addError("Failed reading haplotypes from reference.");
                exit(EXIT_FAILURE);
            }
            readbytes += qin.gcount();
        }
    }

    qrefcurridxglob++;
    qrefcurridxreg++;
}


// read all variants from current stream position until either the position is reached exactly
// (if there are several variants at this position, not all of them might be loaded!)
// or until the closest variant BEFORE that position.
// data is added with "push_back" to referenceFullT.
// returns true, if at least tgt position matches, and only then, the following bools are set accordingly:
// - if ref and alt allele were swapped, refaltswapped is set.
// - if the reverse complements of the alleles match, strandflipped is set.
// - if the tgt ref/alt pair has no counterpart in the reference, refalterror is set.
// qridx at function call is used as starting position to find the variant, we require qridx <= qrefcurridxreg!!
// qridx after return is either the index (relative to region) where the tgt was found or the index of the following variant if not found.
// NOTE: the function may swap entries at positions greater or equal to the returned qridx
//       (happens if tgt was found in a multi-allelic variant, then the matching variant is swapped to qridx, while others are moved to higher indices)
inline bool VCFData::loadAndFindInQref(bcf1_t *tgt, size_t &qridx, bool &refaltswapped, bool &strandflipped, bool &refalterror) {

    // find position in qref (qridx <= qrefcurridxreg!!!)
    while (qridx < Mrefreg && tgt->pos > positionsFullRefRegion[qridx]) {
        if (qridx == qrefcurridxreg) {
            qRefLoadNextVariant();
        }
        // DEBUG
        if (positionsFullRefRegion[qridx] >= 33800342 && positionsFullRefRegion[qridx] <= 33800401) {
            cerr << "ZZZ: Loaded at " << qridx << ": " << positionsFullRefRegion[qridx] << ":" << allelesFullRefRegion[2*qridx] << ":" << allelesFullRefRegion[2*qridx+1] << ":" << (multiAllFlagsFullRefRegion[qridx] ? "ma" : "ba") << endl;
        }
        if (positionsFullRefRegion[qridx] == 34108971) {
            cerr << "XYZ: Loaded at " << qridx << ": " << positionsFullRefRegion[qridx] << ":" << allelesFullRefRegion[2*qridx] << ":" << allelesFullRefRegion[2*qridx+1] << ":" << (multiAllFlagsFullRefRegion[qridx] ? "ma" : "ba") << endl;
        }
        // __DEBUG
        qridx++;
    }
    if (tgt->pos != positionsFullRefRegion[qridx])
        return false; // not found

    // at least positions match

    bool mono = tgt->n_allele == 1; // check, if monomorph
    size_t tryidx = qridx; // continuing now with this variant and try all indices at the same position
    size_t fidx = 0; // if we found the tgt variant, this index is set accordingly
    bcf_unpack(tgt, BCF_UN_STR); // unpack alleles
    string tgtall0(tgt->d.allele[0]);
    string tgtall0rc = reverseComplement(tgt->d.allele[0]);
    string tgtall1;
    string tgtall1rc;
    if (!mono) {
        tgtall1 = string(tgt->d.allele[1]);
        tgtall1rc = reverseComplement(tgt->d.allele[1]);
    }

    // check alleles, also in possible multi-allelic splits
    bool found = false;
    do {
        // load variant (if not yet loaded)
        if (tryidx == qrefcurridxreg)
            qRefLoadNextVariant();
        // DEBUG
        if (positionsFullRefRegion[tryidx] >= 33800342 && positionsFullRefRegion[tryidx] <= 33800401) {
            cerr << "ZZZ: Loaded match at " << tryidx << ": " << positionsFullRefRegion[tryidx] << ":" << allelesFullRefRegion[2*tryidx] << ":" << allelesFullRefRegion[2*tryidx+1] << ":" << (multiAllFlagsFullRefRegion[tryidx] ? "ma" : "ba") << endl;
        }
        if (positionsFullRefRegion[qridx] == 34108971) {
            cerr << "XYZ: Loaded match at " << tryidx << ": " << positionsFullRefRegion[tryidx] << ":" << allelesFullRefRegion[2*tryidx] << ":" << allelesFullRefRegion[2*tryidx+1] << ":" << (multiAllFlagsFullRefRegion[tryidx] ? "ma" : "ba") << endl;
        }
        // __DEBUG

        if (allelesFullRefRegion[2*tryidx].compare(tgtall0) == 0) { // reference alleles match
            if (mono || allelesFullRefRegion[2*tryidx+1].compare(tgtall1) == 0) { // alternative alleles match
                // highest priority match
                refaltswapped = false;
                strandflipped = false;
                found = true;
                fidx = tryidx;
                break; // don't need to try others
            }
        } else if (allowRefAltSwap && allelesFullRefRegion[2*tryidx+1].compare(tgtall0) == 0) { // reference alleles match if swapped
            if (mono || allelesFullRefRegion[2*tryidx].compare(tgtall1) == 0) { // alternative alleles match if swapped
                // second highest priority match:
                // a highest priority match would have exited already, so set without checking, but continue searching
                refaltswapped = true;
                strandflipped = false;
                found = true;
                fidx = tryidx;
            }
        } else if (allowStrandFlip && allelesFullRefRegion[2*tryidx].compare(tgtall0rc) == 0) { // reference alleles match if strands are flipped
            if (mono || allelesFullRefRegion[2*tryidx+1].compare(tgtall1rc) == 0) { // alternative alleles match if strands are flipped
                // third highest priority match:
                // check if we had a second highest match already (highest priority would have exited already)
                if (!(found && refaltswapped && !strandflipped)) {
                    refaltswapped = false;
                    strandflipped = true;
                    found = true;
                    fidx = tryidx;
                }
            }
        } else if (allowRefAltSwap && allowStrandFlip && allelesFullRefRegion[2*tryidx+1].compare(tgtall0rc) == 0) { // reference alleles match if swapped and if strands are flipped
            if (mono || allelesFullRefRegion[2*tryidx].compare(tgtall1rc) == 0) { // alternative alleles match if swapped and if strands are flipped
                // least highest priority match:
                // check if we had a higher priority match before (so, if there was a find at all before)
                if (!found) {
                    refaltswapped = true;
                    strandflipped = true;
                    found = true;
                    fidx = tryidx;
                }
            }
        }
        // no match at all -> try next idx
        tryidx++;
    } while(tryidx < positionsFullRefRegion.size() && tgt->pos == positionsFullRefRegion[tryidx]);

    refalterror = !found;
    if (!found) { // positions match, but no ref/alt pair matches
        // the variant at fidx is not yet loaded
        return true;
    }

    // here, we found the tgt (at fidx), but we probably need to swap multi-allelic splits
    // (since we left the loop above with break, the variant at fidx is loaded!)
    if (fidx != qridx) {
        // swap references at fidx and qridx

        // haplotypes
        // simply need to swap data pointers in underlying BooleanVectors
        auto qrdata = referenceFullT[qridx-currChunkOffset].getData();
        auto fdata = referenceFullT[fidx-currChunkOffset].getData();
        referenceFullT[qridx-currChunkOffset].setData(fdata, hapcapacity, Nrefhaps);
        referenceFullT[fidx-currChunkOffset].setData(qrdata, hapcapacity, Nrefhaps);

        // allele frequencies
        swap(alleleFreqsFullRefRegion[qridx], alleleFreqsFullRefRegion[fidx]);

        // variant positions don't need to be swapped, they are equal!

        // flags for multi-allelic splits -> we swap them here
        // since it may be possible that two variants were located at the same position in the original reference,
        // but not marked as multi-allelic
        bool ma0 = multiAllFlagsFullRefRegion[qridx];
        if (multiAllFlagsFullRefRegion[fidx] != ma0) {
            multiAllFlagsFullRefRegion.set(qridx, multiAllFlagsFullRefRegion[fidx]);
            multiAllFlagsFullRefRegion.set(fidx, ma0);
        }

        // allele strings
        swap(allelesFullRefRegion[2*qridx], allelesFullRefRegion[2*fidx]);
        swap(allelesFullRefRegion[2*qridx+1], allelesFullRefRegion[2*fidx+1]);

        // allele IDs
        swap(variantIDsFullRefRegion[qridx], variantIDsFullRefRegion[fidx]);

        // DEBUG
        if (positionsFullRefRegion[qridx] == 34108971) {
            cerr << "XYZ: Swapped " << qridx << " and " << fidx << ": " << positionsFullRefRegion[qridx] << ":" << allelesFullRefRegion[2*qridx] << ":" << allelesFullRefRegion[2*qridx+1] << ":" << (multiAllFlagsFullRefRegion[qridx] ? "ma" : "ba") << endl;
        }
        // __DEBUG

    }
    return true;
}

// split multi-allelic reference SNPs:
// *ref will be updated to a bi-allelic SNP with the reference and first alternative allele.
// for the other alternative alleles a new record will be created with the same reference allele and put into masplits.
// nSplitted will be incremented by the number of resulting splits (including the original).
// gFilled are the number of haplotypes that had to be filled with reference values in each split.
inline void VCFData::splitSNP(bcf1_t *ref, void **ref_gt, int *mref_gt,
        vector<bcf1_t*> &masplits, size_t &nSplitted, size_t &gFilled) {
    // read genotypes
    int nrefgt = bcf_get_genotypes(ref_hdr, ref, ref_gt, mref_gt);
    // determine allele frequencies
    int nrefs = 0;
    int refnalt = 0;
    vector<int> splitnalts(ref->n_allele-2);
    int *ptr = reinterpret_cast<int*>(*ref_gt);
    for (int g = 0; g < nrefgt; g++, ptr++) {
        if (!bcf_gt_is_missing(*ptr) && *ptr != bcf_int32_vector_end) {
            int alidx = bcf_gt_allele(*ptr);
            switch (alidx) {
            case 0:
                nrefs++;
                break;
            case 1:
                refnalt++;
                break;
            default:
                splitnalts[alidx-2]++;
            }
        }
    }
    int splitnaltsum = 0;
    for (auto splitnalt : splitnalts)
        splitnaltsum += splitnalt;
    if (splitnaltsum == 0) { // in fact, the split is declared multi-allelic, but it is bi-allelic -> don't split
        // but update alleles
        bcf_update_alleles(ref_hdr, ref, (const char**)(ref->d.allele), 2);
        // and formally count this as split
        nSplitted++;
        return;
    }

    // here, we have some multi-allelic content

    // prepare output records: one split for each allele after the second
    for (int i = 2; i < ref->n_allele; i++) {
        masplits.push_back(bcf_dup(ref));
        bcf_unpack(masplits.back(), BCF_UN_STR); // need to unpack the duplicate to be able to edit the alleles
    }
    nSplitted += masplits.size() + 1; // also count the original!

    vector<void*> splitgts(masplits.size()); // space for the genotypes, fill with nullpointers
    vector<int> splitmgts(masplits.size());  // capacity for the genotype's space, fill with zeros

    vector<float> splitafs(masplits.size());
    for (int s = 0; s < (int) masplits.size(); s++) {
        splitafs[s] = ((double) splitnalts[s]) / (nrefs + splitnalts[s]);
    }

    // get genotypes from copies
    for (size_t s = 0; s < masplits.size(); s++)
        bcf_get_genotypes(ref_hdr, masplits[s], &(splitgts[s]), &(splitmgts[s]));

    // reset genotype pointer
    ptr = reinterpret_cast<int*>(*ref_gt);
    // replace each allele index greater 1 with a random value between 0 and 1 in the original, analogue to the splits
    for (int g = 0; g < nrefgt; g++, ptr++) {

        if (!bcf_gt_is_missing(*ptr) && *ptr != bcf_int32_vector_end && bcf_gt_allele(*ptr) != 0) { // some alternative allele
            bool phased = bcf_gt_is_phased(*ptr);
            // change allele in splits
            for (int s = 0; s < (int) masplits.size(); s++) {
                int *splitgt = reinterpret_cast<int*>(splitgts[s]);
                // change in non-corresponding split
                if (s != bcf_gt_allele(*ptr)-2) {
                    int alidx = 0; // choose reference allele
                    splitgt[g] = phased ? bcf_gt_phased(alidx) : bcf_gt_unphased(alidx);
                } else if (bcf_gt_allele(*ptr) > 1) { // translate index to '1'
                    splitgt[g] = phased ? bcf_gt_phased(1) : bcf_gt_unphased(1);
                }
            }
            // change allele in original
            if (bcf_gt_allele(*ptr) > 1) {
                int alidx = 0; // choose reference allele
                *ptr = phased ? bcf_gt_phased(alidx) : bcf_gt_unphased(alidx);
            }

            gFilled += masplits.size();
        }

    }

    // update alleles
    char *tmp[2];
    tmp[0] = ref->d.allele[0]; // reference allele
    int refalidx = 2;
    for (auto split : masplits) {
        tmp[1] = ref->d.allele[refalidx]; // choose right alternative allele
        refalidx++;
        bcf_update_alleles(ref_hdr, split, (const char**) tmp, 2);
        // update ID string to add the alleles
        stringstream ss;
        ss << ref->d.id << "_" << tmp[0] << "_" << tmp[1];
        bcf_update_id(ref_hdr, split, ss.str().c_str()); // change ID name
    }
    bcf_update_alleles(ref_hdr, ref, (const char**) ref->d.allele, 2);
    stringstream ss;
    ss << ref->d.id << "_" << ref->d.allele[0] << "_" << ref->d.allele[1];
    bcf_update_id(ref_hdr, ref, ss.str().c_str()); // change ID name

    // update genotypes
    bcf_update_genotypes(ref_hdr, ref, *ref_gt, nrefgt);
    for (size_t s = 0; s < masplits.size(); s++) {
        bcf_update_genotypes(ref_hdr, masplits[s], splitgts[s], nrefgt);
        free(splitgts[s]); // allocated by HTSlib
    }

}

void VCFData::determineStartPhases(const vector<BooleanVector> &phasedTargets) {
    for (size_t tgt = 0; tgt < Ntarget; tgt++) {
        // for each target find the first heterozygous site at the beginning of the overlap and store the phase
        size_t searchidx = currM-currMOverlap; // beginning of overlap
        bool mat, pat;
        do {
            mat = phasedTargets[2*tgt].at(searchidx);
            pat = phasedTargets[2*tgt+1].at(searchidx);
            searchidx++;
        } while (mat == pat && searchidx < phasedTargets[2*tgt].size());
        // here, we found either the first het site, or we reached the end of the phased target, which means there are no hets in the overlap region
        // anyway, we set the phase to the current mat haplotype
        chunkStartPhase[tgt] = mat;
    }
}

void VCFData::writePhasedConfidences(const vector<fp_type> &totalconfs, const vector<size_t> &ncalls) {
    ofstream ofs(outFileNameConf);
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
    for (size_t nt = DEBUG_TARGET_START; nt <= DEBUG_TARGET_STOP; nt++)
#else
    for (size_t nt = 0; nt < Ntarget; nt++)
#endif
    {
        ofs << targetIDs[nt] << "\t";
        if (haploidsTgt[nt])
            ofs << "haploid" << endl;
        else
            ofs << totalconfs[nt]/ncalls[nt] << endl;
    }
    ofs.close();
}

void VCFData::writeVCFPhased(const vector<BooleanVector> &phasedTargets) {

    string wm = writeMode;
    if (currChunk)
        wm[0] = 'a'; // simply append if writing data after the first chunk
    htsFile *out = hts_open(outFileName.c_str(), wm.c_str());
    if (out == NULL) {
        StatusFile::addError("Could not open file for phasing output.");
        exit(EXIT_FAILURE);
    }

    // set maximum number of helper threads for compression
    hts_set_threads(out, numThreads-1);

    bcf_hdr_t *hdr = tgt_hdr_cpy;
    if (!skipHeader && currChunk == 0) {
        if(bcf_hdr_write(out, hdr)) {
            StatusFile::addError("Failed writing header for phasing output.");
            exit(EXIT_FAILURE);
        }
    }

//    int mtgt_gt = Ntarget*2;
    void *tgt_gt = MyMalloc::malloc(Ntarget * 2 * sizeof(int), "tgt_gt_writePhased"); // mtgt_gt is the size of the tgt_gt space (need void* because of htslib)
    int* tgt_gt_int = reinterpret_cast<int*>(tgt_gt); // update, whenever tgt_gt is updated!!!

    size_t psite = 0; // index of phased sites
    auto rec_it = bcf_pout.begin();
    for (bool isphased : isPhased) {

        if (isphased) { // phased marker

            int64_t bp = (*rec_it)->pos; // zero-based!
            if (startChunkBp <= bp && bp <= endChunkBp) { // check if within output region (without flanks)
                for (size_t i = 0; i < Ntarget; i++) {
                    int *ptr = tgt_gt_int + 2*i;

                    bool missing = targets[i][psite] == Genotype::Miss; // will be slow... sh...
                    // if (missing) only checks if the genotype was missing in the original! Check noImpMissing if the genotype was imputed or if we have to write "missing" again.

                    if (!haploidsTgt[i]) { // diploid sample
                        // We had a lot of filters before, especially filtering multi-allelic markers.
                        // So, for each phased marker '0' represents the reference allele and '1' represents
                        // the alternative allele. We don't need to search for the maximum and minimum index as before...

                        if (missing && noImpMissing) {
                            *ptr = *(ptr+1) = bcf_gt_missing;
                        } else { // not missing or imputed missing -> information is in "phasedTarget"
                            size_t nTargetHap = 2*i;
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                            bool hap0 = false;
                            bool hap1 = false;
                            if (psite < phasedTargets[nTargetHap].size()) { // have we phased this target in DEBUG mode?
                                hap0 = phasedTargets[nTargetHap][psite];
                                hap1 = phasedTargets[nTargetHap+1][psite];
                            }
#else
                            bool hap0 = phasedTargets[nTargetHap][psite];
                            bool hap1 = phasedTargets[nTargetHap+1][psite];
#endif
                            *ptr = bcf_gt_unphased(hap0 ? 1 : 0); // convert allele index to bcf value (phased/unphased does not matter here, but unphased is a bit faster)
                            *(ptr+1) = bcf_gt_phased(hap1 ? 1 : 0); // convert allele index to bcf value (phased)
                        }

                    } else { // haploid sample

                        *(ptr+1) = bcf_int32_vector_end; // only first field carries the information, second is always "vector end"

                        if (missing && noImpMissing) {
                            *ptr = bcf_gt_missing;
                        } else {
                            size_t nTargetHap = 2*i;
                            bool hap = phasedTargets[nTargetHap][psite]; // stored as homozygous diploid anyway
                            *ptr = bcf_gt_unphased(hap ? 1 : 0);
                        }
                    }
                }

                bcf_update_genotypes(hdr, *rec_it, tgt_gt, Ntarget * 2);
                if (bcf_write(out, hdr, *rec_it)) {
                    StatusFile::addError("Failed writing phasing output.");
                    exit(EXIT_FAILURE);
                }

            }
            psite++; // next phased site

        } else { // !isphased: site was not phased; we do not remove phase information anymore
            // basically, we do nothing here, but check if the SNP is in the output region
            int64_t bp = (*rec_it)->pos; // zero-based!
            if (startChunkBp <= bp && bp <= endChunkBp) { // check if within output region
                if (bcf_write(out, hdr, *rec_it)) {
                    StatusFile::addError("Failed writing phasing output.");
                    exit(EXIT_FAILURE);
                }
            }
        }
        bcf_destroy(*rec_it);
        ++rec_it;
    } // for

    MyMalloc::free(tgt_gt);
//    bcf_destroy(rec); // all records are destroyed in the loop
    hts_close(out);
}

// returns the number of sites per file
void VCFData::writeVCFImputedPrepare() {

    // create output files
    bcfouts.clear();
    bcfouts.resize(num_files);
    // first file gets final file name + header, others will be appended to this
    string ofn = currChunk ? string(outFileNameImp).append(".ch").append(to_string(currChunk)) : outFileNameImp;
    bcfouts[0] = hts_open(ofn.c_str(), writeMode.c_str());
    if (bcfouts[0] == NULL) {
        StatusFile::addError("Could not open file for imputation output.");
        exit(EXIT_FAILURE);
    }
    // temporary files
    for (unsigned f = 1; f < num_files; f++) {
        bcfouts[f] = hts_open(string(ofn).append(".").append(to_string(f)).c_str(), writeMode.c_str());
        if (bcfouts[f] == NULL) {
            StatusFile::addError("Could not open temporary file for imputation output.");
            exit(EXIT_FAILURE);
        }
    }

    // set helper threads for compression (only if number of threads > 1)
    if (numThreads > 1) {
        for (unsigned f = 0; f < num_files; f++)
            hts_set_threads(bcfouts[f], num_workers+num_files-1); // equally distribute all available threads to each file, but add at least one.
    }

    // write VCF header
    if (!skipHeader && currChunk == 0) {
        if (bcf_hdr_write(bcfouts[0], imp_hdr)) {
            StatusFile::addError("Failed writing header for imputation output.");
            exit(EXIT_FAILURE);
        }
    } else { // header must be sync'ed although not written
        if (bcf_hdr_sync(imp_hdr) < 0) {
            StatusFile::addError("Failed syncing header for imputation output.");
            exit(EXIT_FAILURE);
        }
    }

    // dosage information (space re-used for each SNP)
    ads_vec.resize(num_workers);
    for (auto &ads : ads_vec)
        ads.resize(Ntarget*2); // for each target and allele the allele dosage (allele dosages are always stored)
    ds_vec.resize(num_workers);
    if (writeGDosage) {
        for (auto &ds : ds_vec)
            ds.resize(Ntarget); // for each target the genotype dosage
    }
    gp_vec.resize(num_workers);
    if (writeProbs) {
        for (auto &gp : gp_vec)
            gp.resize(Ntarget*3); // for each target and genotype the posterior probability
    }

    mtgt_gt_vec.resize(num_workers, Ntarget*2); // mtgt_gt is the size of the tgt_gt space
    tgt_gt_vec.resize(num_workers);
    for (auto &tgt_gt : tgt_gt_vec)
        tgt_gt = MyMalloc::malloc(Ntarget * 2 * sizeof(int), "tgt_gt_writeImputed"); // (need void* because of htslib)

    site_offsets.resize(num_files+1); // field at the end represents the total number of imputation sites (including left overlap!)
    site_offsets[0] = num_sites_per_file[0]; // first file skips the overlap
    for (unsigned f = 1; f <= num_files; f++) {
        site_offsets[f] = site_offsets[f-1] + num_sites_per_file[f]; // keep in mind that the first field in num_sites_.. contains the vars in the overlap
    }
    size_t num_sites_no_ov = site_offsets.back() - num_sites_per_file[0]; // need to remove the first overlap to get the number of sites without overlaps

    // now that we know the exact distribution of sites for each file, we can calculate the optimal bunch size and the number of bunches
    if (bunchsize > divideRounded(num_sites_no_ov, (size_t)num_files)) { // don't exceed the total number of reference variants in current chunk
        bunchsize = divideRounded(num_sites_no_ov, (size_t)num_files);
    }
    nbunches = divideRounded(num_sites_no_ov, bunchsize * num_files);
    // DEBUG
    cerr << "BUNCHES (pre) nb / size: " << nbunches << " / " << bunchsize << endl;
    // optimize bunchsize to have an equal load on all bunches (i.e. reducing the bunchsize to have all bunches approximately the same size)
    bunchsize = divideRounded(num_sites_no_ov, nbunches * num_files);
    // only set to multiple of num_workers if the memory increase would not be too much!
    if (bunchsize > 4*num_workers)
        bunchsize = roundToMultiple(bunchsize, (size_t)num_workers);
    nbunches = divideRounded(num_sites_no_ov, bunchsize * num_files); // don't know if this is required, but it doesn't hurt
    // DEBUG
    cerr << "BUNCHES (post) nb / size: " << nbunches << " / " << bunchsize << endl;

    // create output queues, organized for each file
    size_t qcap = 2 * divideRounded(bunchsize, (size_t)num_workers); // space for two bunches (distributed to each worker queue)
    recqs.resize(num_files);
    for (unsigned f = 0; f < num_files; f++) {
        recqs[f].resize(num_workers);
        for (auto &q : recqs[f])
            q.set_capacity(qcap);
    }

    // start writing threads
    function<void(vector<tbb::concurrent_bounded_queue<bcf1_t*>>&, htsFile*, bcf_hdr_t*, size_t)> wrfunc(std::bind(&VCFData::writeBCFRecords, this, _1, _2, _3, _4));
    for (unsigned f = 0; f < num_files; f++) {
        wrts.emplace_back(wrfunc, std::ref(recqs[f]), std::ref(bcfouts[f]), std::ref(imp_hdr), num_sites_per_file[f+1]);
    }

}


void VCFData::writeVCFImputedBunch(
        unsigned block,
        size_t startidx,
        size_t bsize,
        const vector<BooleanVector> &imputedTargets,
        const vector<vector<float>> &imputedDosages,
        const vector<BooleanVector> &phasedTargets,
        const vector<vector<float>> &phasedDosages) {

    const bcf_hdr_t *bcfhdr = imp_hdr;

    // DEBUG
    atomic<size_t> nfiltered1(0);
    atomic<size_t> nfiltered2(0);
    // DEBUG

    omp_set_num_threads(num_workers);
#pragma omp parallel
    for (size_t bidx = omp_get_thread_num(); bidx < bsize && site_offsets[block] + startidx + bidx < site_offsets[block+1]; bidx += num_workers) { // go through all sites in bunch in current block

        size_t idx = site_offsets[block] + startidx + bidx;
        size_t idxReg = idx + currChunkOffset;

        // DEBUG
        if (idxReg >= 2295058 && idxReg <= 2295064) {
            cerr << "ZZZ: Touched idxReg: " << idxReg << " idx: " << idx << " site_offsets[" << block << "]: " << site_offsets[block] << " startidx: " << startidx << " bidx: " << bidx << endl;
        }
        // __DEBUG

        vector<float> &ads = ads_vec[omp_get_thread_num()];
        vector<float> &ds = ds_vec[omp_get_thread_num()];
        vector<float> &gp = gp_vec[omp_get_thread_num()];

        int* tgt_gt_int = reinterpret_cast<int*>(tgt_gt_vec[omp_get_thread_num()]); // update, whenever tgt_gt_vec[n] is updated!!!

        size_t psite = ptgtIdx[idx]; // index of phased sites (index in phasedTargets)
        size_t isite = bidx;

        bool isimputed = isImputed[idx];

        // check if within output region and not an imputed multi-allelic that the user wants to skip
        int64_t bp = positionsFullRefRegion[idxReg];
        if (startChunkBp <= bp && bp <= endChunkBp && !(isimputed && excludeMultiAllRef && multiAllFlagsFullRefRegion[idxReg]) ) {

            bcf1_t *bcfrec = bcf_init1();

            bcf_float_set_missing(bcfrec->qual);

            // saves time and space if we don't set this
    //        bcf_add_filter(bcfhdr, bcfrec, bcf_hdr_id2int(bcfhdr, BCF_DT_ID, "PASS"));

            bcfrec->rid = 0; // since we only do one chromosome, this has to be the first

            // copy position, variant ID and allele information from the reference
            // NOTE: we only output common variants, so this is completely adequate also for phased and unphased variants
            bcfrec->pos = positionsFullRefRegion[idxReg];
            bcf_update_id(bcfhdr, bcfrec, variantIDsFullRefRegion[idxReg].c_str());
            bcf_update_alleles_str(bcfhdr, bcfrec, string(allelesFullRefRegion[2*idxReg]).append(",").append(allelesFullRefRegion[2*idxReg+1]).c_str());

            float rpaf = alleleFreqsFullRefRegion[idxReg]; // allele frequency in reference panel
            double adsum = 0.0; // sum of allele dosages, rewuired for r2
            double ad2sum = 0.0; // sum of squared allele dosages, required for r2
            int nval = 0; // sum of samples to count for info score calculation
            int ac1 = 0; // allele count for 1-allele
            int an = 0;  // total number of alleles

            if (isimputed) { // site was imputed

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                for (size_t i = DEBUG_TARGET_START; i <= DEBUG_TARGET_STOP; i++)
#else
                for (size_t i = 0; i < Ntarget; i++)
#endif
                {
                    int *ptr = tgt_gt_int + 2*i;
                    if (!haploidsTgt[i]) { // diploid sample

                        size_t nTargetHap = 2*i;
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                        bool hap0 = false;
                        bool hap1 = false;
//                        if (itgtIdx[idx] < phasedTargets[nTargetHap].size()) { // have we imputed this target in DEBUG mode?
                        if (idx < phasedTargets[nTargetHap].size()) { // have we imputed this target in DEBUG mode?
                            hap0 = imputedTargets[nTargetHap][isite];
                            hap1 = imputedTargets[nTargetHap+1][isite];
                        }
#else
                        bool hap0 = imputedTargets[nTargetHap][isite];
                        bool hap1 = imputedTargets[nTargetHap+1][isite];
#endif
                        // dosage
                        float ads0 = imputedDosages[nTargetHap][isite];
                        if (ads0 < 0) { // imputation conflict: set dosage to allele frequency
                            ads0 = rpaf;
                            hap0 = rpaf > 0.5;
                            conflictcnt++;
                        }
                        ads[2*i] = ads0;

                        float ads1 = imputedDosages[nTargetHap+1][isite];
                        if (ads1 < 0) { // imputation conflict: set dosage to allele frequency
                            ads1 = rpaf;
                            hap1 = rpaf > 0.5;
                            conflictcnt++;
                        }
                        ads[2*i+1] = ads1;

                        adsum += ads0 + ads1;
                        ad2sum += ads0*ads0 + ads1*ads1;
                        if (writeGDosage)
                            ds[i] = ads0 + ads1;
                        if (writeProbs) {
                            double gp0 = (1.0 - ads0) * (1.0 - ads1);
                            double gp2 = ads0 * ads1;
                            double gp1 = 1.0 - gp0 - gp2;
                            gp[3*i]   = gp0;
                            gp[3*i+1] = gp1;
                            gp[3*i+2] = gp2;
                        }

                        // allele count
                        if (hap0)
                            ac1++;
                        if (hap1)
                            ac1++;
                        an += 2;
                        nval += 2;

                        // genotype
                        *ptr = bcf_gt_unphased(hap0 ? 1 : 0); // convert allele index to bcf value (phased/unphased does not matter here, but unphased is a bit faster)
                        *(ptr+1) = bcf_gt_phased(hap1 ? 1 : 0); // convert allele index to bcf value (phased)

                    } else { // haploid sample

                        size_t nTargetHap = 2*i;
                        bool hap = imputedTargets[nTargetHap][isite]; // stored as homozygous diploid anyway

                        // dosage
                        float ads0 = imputedDosages[nTargetHap][isite];
                        if (ads0 < 0) { // imputation conflict: set dosage to allele frequency
                            ads0 = rpaf;
                            hap = rpaf > 0.5;
                            conflictcnt++;
                        }
                        ads[2*i] = ads0;
                        bcf_float_set(&ads[2*i+1], bcf_float_vector_end);

                        adsum += ads0;
                        ad2sum += ads0*ads0;
                        if (writeGDosage)
                            ds[i] = ads0;
                        if (writeProbs) {
                            double gp0 = (1.0 - ads0);
                            double gp2 = ads0;
                            gp[3*i] = gp0;
                            gp[3*i+1] = gp2;
                            bcf_float_set(&gp[3*i+2], bcf_float_vector_end);
                        }

                        // allele count
                        if (hap)
                            ac1++;
                        an++;
                        nval++;

                        // genotype
                        *ptr = bcf_gt_unphased(hap ? 1 : 0);
                        *(ptr+1) = bcf_int32_vector_end; // only first field carries the information, second is always "vector end"

                    }
                } // END for all targets

            } else { // not imputed

#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                for (size_t i = DEBUG_TARGET_START; i <= DEBUG_TARGET_STOP; i++)
#else
                for (size_t i = 0; i < Ntarget; i++)
#endif
                {
                    int *ptr = tgt_gt_int + 2*i;

                    bool missing = targets[i][psite] == Genotype::Miss; // will be slow... sh...
                    // if (missing) only checks if the genotype was missing in the original! Check noImpMissing if the genotype was imputed or if we have to write "missing" again.

                    if (!haploidsTgt[i]) { // diploid sample
                        // We had a lot of filters before, especially filtering multi-allelic markers.
                        // So, for each phased marker '0' represents the reference allele and '1' represents
                        // the alternative allele. We don't need to search for the maximum and minimum index as before...

                        size_t nTargetHap = 2*i;
                        bool hap0, hap1;
                        float ads0, ads1;
                        if (missing || overwriteCalls) { // site was imputed or we take the imputation result anyway
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                            hap0 = false;
                            hap1 = false;
//                                if (itgtIdx[idx] < phasedTargets[nTargetHap].size()) { // have we imputed this target in DEBUG mode?
                            if (idx < phasedTargets[nTargetHap].size()) { // have we imputed this target in DEBUG mode?
                                hap0 = imputedTargets[nTargetHap][isite];
                                hap1 = imputedTargets[nTargetHap+1][isite];
                            }
#else
                            hap0 = imputedTargets[nTargetHap][isite];
                            hap1 = imputedTargets[nTargetHap+1][isite];
#endif
                            ads0 = imputedDosages[nTargetHap][isite];
                            if (ads0 < 0) { // imputation conflict: set dosage to allele frequency
                                ads0 = rpaf;
                                hap0 = rpaf > 0.5;
                                conflictcnt++;
                            }
                            ads[2*i] = ads0;

                            ads1 = imputedDosages[nTargetHap+1][isite];
                            if (ads1 < 0) { // imputation conflict: set dosage to allele frequency
                                ads1 = rpaf;
                                hap1 = rpaf > 0.5;
                                conflictcnt++;
                            }
                            ads[2*i+1] = ads1;

                        } else { // not missing -> information is in "phasedTarget"
#if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                            hap0 = false;
                            hap1 = false;
                            if (psite < phasedTargets[nTargetHap].size()) { // have we phased this target in DEBUG mode?
                                hap0 = phasedTargets[nTargetHap][psite];
                                hap1 = phasedTargets[nTargetHap+1][psite];
                            }
#else
                            hap0 = phasedTargets[nTargetHap][psite];
                            hap1 = phasedTargets[nTargetHap+1][psite];
#endif
                            ads0 = phasedDosages[nTargetHap][psite];
                            ads1 = phasedDosages[nTargetHap+1][psite];

                            // take the imputation result if this would improve the dosage
                            if (improveCalls) {
                                float idos0 = imputedDosages[nTargetHap][isite];
                                float idos1 = imputedDosages[nTargetHap+1][isite];
                                float pdos0 = ads0;
                                float pdos1 = ads1;
                                // take allele frequency on imputation conflict
                                if (idos0 < 0)
                                    idos0 = rpaf;
                                if (idos1 < 0)
                                    idos1 = rpaf;
                                // "normalize" -> the closer to zero the better
                                if (idos0 > 0.5)
                                    idos0 = 1.0 - idos0;
                                if (idos1 > 0.5)
                                    idos1 = 1.0 - idos1;
                                if (pdos0 > 0.5)
                                    pdos0 = 1.0 - pdos0;
                                if (pdos1 > 0.5)
                                    pdos1 = 1.0 - pdos1;
                                float improve0 = pdos0 - idos0;
                                float improve1 = pdos1 - idos1;

                                if (improve0+improve1 > 0) {
                                    // some statistics
                                    if (imputedDosages[nTargetHap][isite] < 0)
                                        conflictcnt++;
                                    if (imputedDosages[nTargetHap+1][isite] < 0)
                                        conflictcnt++;
                                    if (hap0 != imputedTargets[nTargetHap][isite] && hap1 != imputedTargets[nTargetHap+1][isite]) {
                                        if (hap0 != hap1) {
                                            improvecnthet2het++;
                                            totalimprovehet2het += improve0+improve1;
                                        } else {
                                            improvecnthom2hom++;
                                            totalimprovehom2hom += improve0+improve1;
                                        }
                                    } else if (hap0 != imputedTargets[nTargetHap][isite] || hap1 != imputedTargets[nTargetHap+1][isite]) {
                                        if (hap0 == hap1) {
                                            improvecnthom2het++;
                                            totalimprovehom2het += improve0+improve1;
                                        } else {
                                            improvecnthet2hom++;
                                            totalimprovehet2hom += improve0+improve1;
                                        }
                                    } else {
                                        improvecntnochange++;
                                        totalimprovenochange += improve0+improve1;
                                    }

                                    // change the alleles and dosages
                                    hap0 = imputedTargets[nTargetHap][isite];
                                    hap1 = imputedTargets[nTargetHap+1][isite];
                                    ads0 = hap0 ? 1.0 - idos0 : idos0;
                                    ads1 = hap1 ? 1.0 - idos1 : idos1;
                                }
                            }

                            ads[2*i] = ads0;
                            ads[2*i+1] = ads1;
                        }
                        // dosage
                        adsum += ads0 + ads1;
                        ad2sum += ads0*ads0 + ads1*ads1;
                        if (writeGDosage)
                            ds[i] = ads0 + ads1;
                        if (writeProbs) {
                            double gp0 = (1.0 - ads0) * (1.0 - ads1);
                            double gp2 = ads0 * ads1;
                            double gp1 = 1.0 - gp0 - gp2;
                            gp[3*i]   = gp0;
                            gp[3*i+1] = gp1;
                            gp[3*i+2] = gp2;
                        }
                        // allele count
                        if (hap0)
                            ac1++;
                        if (hap1)
                            ac1++;
                        an += 2;
                        nval += 2;
                        // genotype
                        *ptr = bcf_gt_unphased(hap0 ? 1 : 0); // convert allele index to bcf value (phased/unphased does not matter here, but unphased is a bit faster)
                        *(ptr+1) = bcf_gt_phased(hap1 ? 1 : 0); // convert allele index to bcf value (phased)

                    } else { // haploid sample

                        size_t nTargetHap = 2*i;
                        bool hap;
                        if (missing || overwriteCalls) { // site was imputed
                            float ads0;
                            hap = imputedTargets[nTargetHap][isite]; // stored as homozygous diploid anyway
                            ads0 = imputedDosages[nTargetHap][isite];
                            if (ads0 < 0) { // imputation conflict: set dosage to allele frequency
                                ads0 = rpaf;
                                hap = rpaf > 0.5;
                                conflictcnt++;
                            }
                            bcf_float_set(&ads[2*i+1], bcf_float_vector_end);
                            adsum += ads0;
                            ad2sum += ads0*ads0;
                            if (writeGDosage)
                                ds[i] = ads0;
                            if (writeProbs) {
                                double gp0 = (1.0 - ads0);
                                double gp2 = ads0;
                                gp[3*i] = gp0;
                                gp[3*i+1] = gp2;
                                bcf_float_set(&gp[3*i+2], bcf_float_vector_end);
                            }
                            // allele count
                            if (hap)
                                ac1++;
                            an++;
                            nval++;
                            ads[2*i] = ads0;
                        } else {
                            hap = phasedTargets[nTargetHap][psite]; // stored as homozygous diploid anyway
                            bcf_float_set(&ads[2*i], bcf_float_missing);
                            bcf_float_set(&ads[2*i+1], bcf_float_vector_end);
                        }

                        // genotype
                        *ptr = bcf_gt_unphased(hap ? 1 : 0);
                        *(ptr+1) = bcf_int32_vector_end; // only first field carries the information, second is always "vector end"
                    } // END haploid sample
                } // END for all targets

            } // END !isImputed

            // calculate r2 score from imputation
            // the default will be 1.0 for TYPED_ONLY SNPs, all others have generated dosage values which will be used for score calculation
            // equations as in minimac4, calculation optimized
            double r2 = 1.0;
            double af = 0.0;
            double evar = 0.0;
            double ovar = 0.0;
            if (nval) { // values for r2 calculation exist
                af = adsum / nval;
                evar = af * (1.0 - af);
                ovar = (ad2sum - adsum * af) / nval;
                if (evar > 0.0 && ovar > 0.0)
                    r2 = ovar / evar;
                else
                    r2 = 0.0;
            }

            // convert to float
            // (needed for writing vcf/bcf tag)
            // (filtering is also done according to the float value to be consistent with the output)
            float r2f = r2;
            float aff = af;
            float maff = aff > 0.5 ? 1-aff : aff;

            // only proceed if score exceeds filter threshold(s)
            if (r2f < impR2filter || maff < impMAFfilter) {
                // add null pointer to output queue to indicate that this record was filtered
                recqs[block][(startidx+bidx) % num_workers].push(NULL);
                bcf_destroy(bcfrec); // not required
            } else {
                // update and enqueue record for writing

                bcf_update_info_float(bcfhdr, bcfrec, "RefPanelAF", &rpaf, 1); // allele frequency in reference panel
                bcf_update_info_float(bcfhdr, bcfrec, "AF", &aff, 1); // estimated allele frequency from allele dosages
                bcf_update_info_float(bcfhdr, bcfrec, "MAF", &maff, 1); // estimated minor allele frequency from allele dosages
                bcf_update_info_float(bcfhdr, bcfrec, "R2", &r2f, 1); // minimac4 imputation info score
                bcf_update_info_int32(bcfhdr, bcfrec, "AC", &ac1, 1); // allele count
                bcf_update_info_int32(bcfhdr, bcfrec, "AN", &an, 1);  // total alleles
                if (!isimputed) { // set flag if typed
                    bcf_update_info_flag(bcfhdr, bcfrec, "TYPED", NULL, 1); // phased and/or imputed anyway
                }

    #if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                bcf_update_genotypes(bcfhdr, bcfrec, tgt_gt_vec[omp_get_thread_num()], 2*(DEBUG_TARGET_STOP-DEBUG_TARGET_START+1));
    #else
                bcf_update_genotypes(bcfhdr, bcfrec, tgt_gt_vec[omp_get_thread_num()], 2*Ntarget);
    #endif

    #if defined DEBUG_TARGET || defined DEBUG_TARGET_LIGHT || defined DEBUG_TARGET_SILENT
                if (writeADosage)
                    bcf_update_format_float(bcfhdr, bcfrec, "ADS", ads.data(), 2*(DEBUG_TARGET_STOP-DEBUG_TARGET_START+1));
                if (writeGDosage)
                    bcf_update_format_float(bcfhdr, bcfrec, "DS", ds.data(), DEBUG_TARGET_STOP-DEBUG_TARGET_START+1);
                if (writeProbs)
                    bcf_update_format_float(bcfhdr, bcfrec, "GP", gp.data(), 3*(DEBUG_TARGET_STOP-DEBUG_TARGET_START+1));
    #else
                if (writeADosage)
                    bcf_update_format_float(bcfhdr, bcfrec, "ADS", ads.data(), 2*Ntarget);
                if (writeGDosage)
                    bcf_update_format_float(bcfhdr, bcfrec, "DS", ds.data(), Ntarget);
                if (writeProbs)
                    bcf_update_format_float(bcfhdr, bcfrec, "GP", gp.data(), 3*Ntarget);
    #endif
    //            }

                // add record to writing queue
                recqs[block][(startidx+bidx) % num_workers].push(bcfrec); // ensure correct order when writing data!

            } // END filter

        } else { // not in output region or skipped multi-allelic

            recqs[block][(startidx+bidx) % num_workers].push(NULL);

            // DEBUG
            if (bp < startChunkBp) nfiltered1++;
            if (bp > endChunkBp) nfiltered2++;
            // __DEBUG

        } // END if not in output region

    } // END through all sites from bunch

    // DEBUG
    if (nfiltered1 || nfiltered2)
        cerr << "YYY: block: " << block << " startidx: " << startidx << " nfiltered1: " << nfiltered1 << " nfiltered2: " << nfiltered2 << endl;
}

void VCFData::writeVCFImputedClose() {
    // wait for writing threads to finish
    for (auto &wrt : wrts)
        wrt.join();
    wrts.clear();

    for (auto &tgt_gt : tgt_gt_vec)
        MyMalloc::free(tgt_gt);

    for (auto &bcfout : bcfouts)
        hts_close(bcfout);
}

void VCFData::writeVCFImputedConcat() {
    // concatenate temporary files
    cout << "Concatenating temporary files..." << endl;
    string ofn = currChunk ? string(outFileNameImp).append(".ch").append(to_string(currChunk)) : outFileNameImp;
    vector<string> tmpfiles;
    tmpfiles.push_back(ofn);
    for (unsigned f = 1; f < num_files; f++)
        tmpfiles.push_back(string(ofn).append(".").append(to_string(f)));
    concatFiles(tmpfiles);
}

void VCFData::combineChunks() const {
    if (nChunks > 1) {
        cout << "Combining chunk files..." << endl;
//        StatusFile::updateStatus(0, "Combine chunks");
        vector<string> chunkfiles;
        chunkfiles.push_back(outFileNameImp);
        for (int chunk = 1; chunk < nChunks; chunk++)
            chunkfiles.push_back(string(outFileNameImp).append(".ch").append(to_string(chunk)));
        concatFiles(chunkfiles);
    }
}

inline void VCFData::appendVersionToBCFHeader(bcf_hdr_t *hdr) const {
    stringstream version;
    version << "##EagleImpVersion=" << prog_version << " htslib=" << hts_version() << endl;
    bcf_hdr_append(hdr, version.str().c_str());

    stringstream command;
    command << "##EagleImpCommand=" << command_argv[0];
    for (int i = 1; i < command_argc; i++)
        command << " " << command_argv[i];
    command << endl;
    bcf_hdr_append(hdr, command.str().c_str());
}

inline void VCFData::concatFiles(const vector<string>& filenames) const {
//    ofstream dest(filenames[0].c_str(), ios_base::binary | ios_base::app); // open for appending
//    for (unsigned f = 1; f < filenames.size(); f++) {
//        ifstream src(filenames[f].c_str(), ios_base::binary);
//        dest << src.rdbuf();
//        src.close(); // append to first file
//        remove(filenames[f].c_str()); // delete temporary file
//    }
//    dest.close();
}

inline string VCFData::getOutputSuffix() {
    if (!writeMode.compare("wbu") || !writeMode.compare("wb"))
        return string(".bcf");
    if (!writeMode.compare("wz"))
        return string(".vcf.gz");
    return string(".vcf"); // w
}

inline size_t VCFData::getNumVariantsFromIndex(const string &vcffilename) {
    size_t numvars = 0;

    htsFile* file = hts_open(vcffilename.c_str(), "r");
    if (!file) {
        StatusFile::addError("Could not read VCF file.");
        exit(EXIT_FAILURE);
    }
    bcf_hdr_t *hdr = bcf_hdr_read(file);
    if (!hdr) {
        StatusFile::addError("Could not read header from VCF file.");
        exit(EXIT_FAILURE);
    }

    tbx_t *tbx = NULL;
    hts_idx_t *idx = NULL;

    if (hts_get_format(file)->format == vcf) {
        tbx = tbx_index_load(vcffilename.c_str());
        if (!tbx) {
            StatusFile::addError("Could not read index file for VCF reference.");
            exit(EXIT_FAILURE);
        }
    }
    else if (hts_get_format(file)->format == bcf)
    {
        idx = bcf_index_load(vcffilename.c_str());
        if (!idx) {
            StatusFile::addError("Could not read index file for VCF reference.");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        StatusFile::addError("Could not detect the reference file type as VCF or BCF.");
        exit(EXIT_FAILURE);
    }

    int nseq;
    // just need this to get the number of sequences, so we can cleanup afterwards immediately
    const char** seq = tbx ? tbx_seqnames(tbx, &nseq) : bcf_index_seqnames(idx, hdr, &nseq);
    free(seq); // allocated by HTSlib
    for (int i = 0; i < nseq; i++) {
        uint64_t nrecords, unmapped;
        hts_idx_get_stat(tbx ? tbx->idx : idx, i, &nrecords, &unmapped);
        numvars += nrecords;
    }

    hts_close(file);
    bcf_hdr_destroy(hdr);
    if (tbx)
        tbx_destroy(tbx);
    if (idx)
        hts_idx_destroy(idx);

    return numvars;
}

/*static*/ int VCFData::chrNameToNumber(const string &name, string &chrlit, bool *literally) {
    if (literally) *literally = false;
	chrlit.assign(name); // simply copy the literal name
    size_t start = 0;
    if (!name.substr(0,3).compare("chr")) {
        // if you named your chromosome "chr"... well, your choice... else, we assume that you prefixed your chromosome name with "chr".
        if (name.size() > 3) {
			start = 3;
        } // else: the name of the chromosome is "chr"
    }
    // we map X/Y to 23,24; other literal forms are kept and mapped to chromosome number 0
    if (start < name.size() && (name[start] < '0' || name[start] > '9')) { // not a number
        if (literally) *literally = true;
        if (start+1 == name.size() && (name[start] == 'X' || name[start] == 'x' || name[start] == 'Y' || name[start] == 'y'))
        	return name[start] == 'X' || name[start] == 'x' ? CHRX : CHRY;
        else
        	return 0; // neither X nor Y, so some other literal. Internally, this is handled with the chromosome number 0.
    }
    // if the name starts with a digit, we map the name to the represented number. if literals follow, we set the literally flag again
    int tmp = 0;
    if (start < name.size() && name[start] >= '0' && name[start] <= '9') {
    	size_t litstart;
        tmp = stoi(name.substr(start), &litstart);
        if (litstart+start < name.size()) { // there is something following the number
        	if (literally) *literally = true;
        }
    // we allow any kind of chromosome now
//    if (tmp < 1 || tmp > 24) {
//        StatusFile::addError("Invalid chromosome! Valid chromosomes are 1-24 or X,Y, preliminary \"chr\" is allowed.");
//        exit(EXIT_FAILURE);
//    }
    }
    return tmp;
}

void VCFData::writeBCFRecords(vector<tbb::concurrent_bounded_queue<bcf1_t*>> &recq, htsFile *bcfout, bcf_hdr_t *bcfhdr, size_t numrecs) {
    ThreadUtils::setThreadName("eglimp write");
//    // DEBUG
//    int maxqsize = 0;
//    // __DEBUG
    auto currq = recq.begin();
    for (size_t cnt = 0; cnt < numrecs; cnt++) {

        bcf1_t *bcfrec;
        currq->pop(bcfrec);
        if (bcfrec) { // do not process null pointers! (these were filtered records)
            if (bcf_write(bcfout, bcfhdr, bcfrec)) {
                StatusFile::addError("Failed writing imputation output.");
                exit(EXIT_FAILURE);
            }
            bcf_destroy(bcfrec);
        }
        currq++;
        if (currq == recq.end())
            currq = recq.begin();
//        // DEBUG
//        int qsize = currq->size();
//        if (qsize > maxqsize)
//            maxqsize = qsize;
    }
//    // DEBUG
//    cout << "writeBCFRecords(): Maximum queue size: " << maxqsize << endl;
}

void VCFData::printSummary() const {
    cout << "\n----------------------------------------------------------------\n--- ";
    stringstream s;
	s << "<h4>Summary</h4>";
	StatusFile::setContext(s.str());
	cout << "Summary:" << endl;
    cout << "----------------------------------------------------------------" << endl;

    StatusFile::addInfo("<p class='pinfo'><b>" + to_string(globalstats.M) + " variants in both target and reference were used for phasing.</b></p>", false);
    size_t Mreftmp = loadQuickRef ? globalstats.Mref-globalstats.MrefMultiAllreg : globalstats.Mref;
    if (doImputation) {
        StatusFile::addInfo("<p class='pinfo'>" + to_string(Mreftmp-globalstats.M) + " variants exclusively in reference were imputed.<br>", false);
        StatusFile::addInfo("<b>Imputation output contains " +to_string(Mreftmp) + " variants.</b></p>", false);
//        // DEBUG
//        cout << "M: " << globalstats.M << endl;
//        cout << "Mref: " << globalstats.Mref << endl;
//        cout << "MrefMultiAllreg: " << globalstats.MrefMultiAllreg << endl;
//        cout << "MmultiAllelicRefTgt: " << globalstats.MmultiAllelicRefTgt << endl;
//        cout << "MmultiAllelicRefOnly: " << globalstats.MmultiAllelicRefOnly << endl;
//        // __DEBUG

    }

    // for YAML messages, placed only if yaml is set
    stringstream yamlinfo;
    yamlinfo << "\n";
    yamlinfo << "    Common variants: " << globalstats.M << "\n";
    if (doImputation) {
        yamlinfo << "    Exclusive reference variants: " << (Mreftmp - globalstats.M) << "\n";
        yamlinfo << "    Imputation output variants: " << Mreftmp << "\n";
    }

    // determine SNP rate (only if we do phasing)
    if (!skipPhasing && nHaploidsTgt != Ntarget) {
        size_t physRange = chrBpsReg.back()-chrBpsReg[0];
        fp_type cMrange = cMs.back() - cM0;

        StatusFile::addInfo("<table>", false);
        StatusFile::addInfo("<tr><td>Physical distance range:        </td><td>" + to_string(physRange) + " base pairs</td></tr>", false);
        StatusFile::addInfo("<tr><td>Genetic distance range:         </td><td>" + to_string(cMrange) + " cM</td></tr>", false);
        StatusFile::addInfo("<tr><td>Average #SNPs per cM in target: </td><td>" + to_string((int)(globalstats.M/cMrange + 0.5)) + "</td></tr></table>", false);

        yamlinfo << "    Physical distance range: " << physRange << " bp\n";
        yamlinfo << "    Genetic distance range: " << cMrange << " cM\n";
        yamlinfo << "    Average SNPs per cM in target: " << (int)(globalstats.M/cMrange + 0.5) << "\n";
    }

    StatusFile::addInfo("<p class='pinfo'>", false);
    if (globalstats.numRefAltSwapAndStrandFlip) {
        StatusFile::addInfo("  REF/ALT were swapped AND strands were flipped in " + to_string(globalstats.numRefAltSwapAndStrandFlip) + " variants.<br>", false);
        yamlinfo << "    REF ALT swap strand flips: " << globalstats.numRefAltSwapAndStrandFlip << "\n";
    }
    if (globalstats.numRefAltSwaps) {
        StatusFile::addInfo("  REF/ALT were swapped in " + to_string(globalstats.numRefAltSwaps) + " variants.<br>", false);
        yamlinfo << "    REF ALT swaps: " << globalstats.numRefAltSwaps << "\n";
    }
    if (globalstats.numStrandFlips) {
        StatusFile::addInfo("  Strands were flipped in " + to_string(globalstats.numStrandFlips) + " variants.<br>", false);
        yamlinfo << "    Strand flips: " << globalstats.numStrandFlips << "\n";
    }

    StatusFile::addInfo("  Dropped " + to_string(globalstats.MtargetOnly) + " variants not found in reference.<br>", false);
    yamlinfo << "    Dropped target only variants: " << globalstats.MtargetOnly << "\n";
    if (globalstats.MmultiAllelicTgt) {
        StatusFile::addInfo("  Dropped " + to_string(globalstats.MmultiAllelicTgt) + " multi-allelic variants in target.<br>", false);
        yamlinfo << "    Dropped multi-allelic target variants: " << globalstats.MmultiAllelicTgt << "\n";
    }
    if (globalstats.MmultiAllelicRefTgt) {
        if (excludeMultiAllRef) {
            StatusFile::addInfo("  Dropped " + to_string(globalstats.MmultiAllelicRefTgt) + " variants bi-allelic in target but multi-allelic in reference.<br>", false);
            yamlinfo << "    Dropped multi-allelic reference variants: " << globalstats.MmultiAllelicRefTgt << "\n";
        } else {
            if (!loadQuickRef){
                StatusFile::addInfo("  Split " + to_string(globalstats.MmultiAllelicRefTgt) + " variants bi-allelic in target but multi-allelic in reference into " + to_string(globalstats.MmultiSplittedTgt)
                        + " bi-allelic reference variants,<br>\n"
                        + "  "  + to_string(globalstats.GmultiFilledTgt) + " haplotypes filled with reference alleles.<br>", false);
                yamlinfo << "    Multi-allelic splits in common variants:\n";
                yamlinfo << "      Common multi-allelic reference variants before split: " << globalstats.MmultiAllelicRefTgt << "\n";
                yamlinfo << "      Resulting bi-allelic reference variants: " << globalstats.MmultiSplittedTgt << "\n";
                yamlinfo << "      Haplotypes filled with reference alleles: " << globalstats.GmultiFilledTgt << "\n";
            } else {
                StatusFile::addInfo("  " + to_string(globalstats.MmultiAllelicRefTgt) + " variants bi-allelic in target were multi-allelic in reference.<br>", false);
                yamlinfo << "    Common reference variants from multi-allelic splits: " << globalstats.MmultiAllelicRefTgt << "\n";
            }
        }
    }
    if (globalstats.MmonomorphicRefTgt) {
        StatusFile::addInfo("  Dropped " + to_string(globalstats.MmonomorphicRefTgt) + " variants bi-allelic in target but monomorphic in reference.<br>", false);
        yamlinfo << "    Dropped target variants monomorphic in reference: " << globalstats.MmonomorphicRefTgt << "\n";
    }
    if (globalstats.MrefAltError) {
        StatusFile::addInfo("  Dropped " + to_string(globalstats.MrefAltError) + " variants with allele mismatches.<br>", false);
        yamlinfo << "    Dropped allele mismatched variants: " << globalstats.MrefAltError << "\n";
    }
    if (useExclude) {
        StatusFile::addInfo("  Dropped " + to_string(globalstats.Mexclude) + " variants from target based on --vcfExclude.<br>", false);
        yamlinfo << "    User excluded variants: " << globalstats.Mexclude << "\n";
    }

    if (globalstats.Munphased) {
        StatusFile::addInfo("  Total unphased variants in phasing output: " + to_string(globalstats.Munphased) + "<br>", false);
        yamlinfo << "    Unphased variants after phasing: " << globalstats.Munphased << "\n";
    }

    StatusFile::addInfo("</p><p class='pinfo'>", false);
    size_t mrefonly = loadQuickRef ? (globalstats.Mref-globalstats.M) : globalstats.MrefOnly;
    StatusFile::addInfo("  " + to_string(mrefonly) + " variants in reference but not in target.<br>", false);
    yamlinfo << "    Reference-only variants: " << mrefonly << "\n";
    if (doImputation) {
        if (excludeMultiAllRef) {
            StatusFile::addInfo("  Excluded " + to_string(globalstats.MmultiAllelicRefOnly) + " reference-only multi-allelic variants from imputation.<br>", false);
            yamlinfo << "    Excluded reference-only multi-allelic variants: " << globalstats.MmultiAllelicRefOnly << "\n";
        } else if (globalstats.MmultiAllelicRefOnly) {
//            if (!loadQuickRef){
//                StatusFile::addInfo("  Split " + to_string(globalstats.MmultiAllelicRefOnly) + " reference-only multi-allelic variants into "
//                        + to_string(globalstats.MmultiSplittedRefOnly) + " bi-allelic reference variants,<br>\n"
//                        + "    " + to_string(globalstats.GmultiFilledRefOnly) + " haplotypes filled with reference alleles.<br>", false);
//                yamlinfo << "    Multi-allelic splits in reference-only variants:\n";
//                yamlinfo << "      Reference-only multi-allelic variants before split: " << globalstats.MmultiAllelicRefOnly << "\n";
//                yamlinfo << "      Resulting bi-allelic reference variants: " << globalstats.MmultiSplittedRefOnly << "\n";
//                yamlinfo << "      Haplotypes filled with reference alleles: " << globalstats.GmultiFilledRefOnly << "\n";
//            } else {
            StatusFile::addInfo("  " + to_string(globalstats.MmultiAllelicRefOnly) + " bi-allelic variants in reference only that originate from multi-allelic splits.<br>", false);
            yamlinfo << "    Reference-only bi-allelic variants from multi-allelic splits: " << globalstats.MmultiAllelicRefOnly << "\n";
//            }
        } else {
            StatusFile::addInfo("  No multi-allelic reference-only variants.<br>", false);
            // we skip this output in YAML
        }
    }
    StatusFile::addInfo("</p><p class='pinfo'>", false);

    if (globalstats.MwithMissingRef) {
        StatusFile::addInfo("  Missing genotypes in phasing reference: " + to_string(globalstats.GmissingRef) + " in " + to_string(globalstats.MwithMissingRef) + " variants.<br>", false);
        yamlinfo << "    Missing genotypes in phasing reference:\n";
        yamlinfo << "      Number: " << globalstats.GmissingRef << "\n";
        yamlinfo << "      Fraction: " << ((globalstats.GmissingRef / (double) globalstats.M) / Nref) << "\n";
        yamlinfo << "      Affected variants: " << globalstats.MwithMissingRef << "\n";
        yamlinfo << "      Affected variants fraction: " << (globalstats.MwithMissingRef / (double) globalstats.M) << "\n";
        StatusFile::addWarning("Reference for phasing contains missing genotypes (set randomly according to allele frequency).<br>\n"
                "   Variants with missing data: " + to_string(globalstats.MwithMissingRef)
                + " Fraction: " + to_string(globalstats.MwithMissingRef / (double) globalstats.M) + "<br>\n"
                "   Missing genotypes: " + to_string(globalstats.GmissingRef)
                + " Fraction: " + to_string((globalstats.GmissingRef / (double) globalstats.M) / Nref));
    }
    if (globalstats.MwithUnphasedRef) {
        StatusFile::addInfo("  Unphased genotypes in phasing reference: " + to_string(globalstats.GunphasedRef) + " in " + to_string(globalstats.MwithUnphasedRef) + " variants.<br>", false);
        yamlinfo << "    Unphased genotypes in phasing reference:\n";
        yamlinfo << "      Number: " << globalstats.GunphasedRef << "\n";
        yamlinfo << "      Fraction: " << ((globalstats.GunphasedRef / (double) globalstats.M) / Nref) << "\n";
        yamlinfo << "      Affected variants: " << globalstats.MwithUnphasedRef << "\n";
        yamlinfo << "      Affected variants fraction: " << (globalstats.MwithUnphasedRef / (double) globalstats.M) << "\n";
        StatusFile::addWarning("Reference for phasing contains unphased genotypes (set to random phase).<br>\n"
                "   Variants with unphased data: " + to_string(globalstats.MwithUnphasedRef)
                + " Fraction: " + to_string(globalstats.MwithUnphasedRef / (double) globalstats.M) + "<br>\n"
                "   Unphased genotypes: " + to_string(globalstats.GunphasedRef)
                + " Fraction: " + to_string((globalstats.GunphasedRef / (double) globalstats.M) / Nref));
    }

    if (doImputation) {
        if (globalstats.MmissingRefOnly) {
            StatusFile::addInfo("  Missing genotypes in imputation reference: " + to_string(globalstats.GmissingRefOnly) + " in " + to_string(globalstats.MmissingRefOnly) + " variants.<br>", false);
            yamlinfo << "    Missing genotypes in imputation reference:\n";
            yamlinfo << "      Number: " << globalstats.GmissingRefOnly << "\n";
            yamlinfo << "      Fraction: " << ((globalstats.GmissingRefOnly / (double) globalstats.Mref) / Nref) << "\n";
            yamlinfo << "      Affected variants: " << globalstats.MmissingRefOnly << "\n";
            yamlinfo << "      Affected variants fraction: " << (globalstats.MmissingRefOnly / (double) globalstats.Mref) << "\n";
            StatusFile::addWarning("Reference for imputation contains missing genotypes (set randomly according to allele frequency).<br>\n"
                    "   Variants with missing data: " + to_string(globalstats.MmissingRefOnly)
                    + " Fraction: " + to_string(globalstats.MmissingRefOnly / (double) globalstats.Mref) + "<br>\n"
                    "   Missing genotypes: " + to_string(globalstats.GmissingRefOnly)
                    + " Fraction: " + to_string((globalstats.GmissingRefOnly / (double) globalstats.Mref) / Nref));
        }
        if (globalstats.MunphasedRefOnly) {
            StatusFile::addInfo("  Unphased genotypes in imputation reference: " + to_string(globalstats.GunphasedRefOnly) + " in " + to_string(globalstats.MunphasedRefOnly) + " variants.<br>", false);
            yamlinfo << "    Unphased genotypes in imputation reference:\n";
            yamlinfo << "      Number: " << globalstats.GunphasedRefOnly << "\n";
            yamlinfo << "      Fraction: " << ((globalstats.GunphasedRefOnly / (double) globalstats.Mref) / Nref) << "\n";
            yamlinfo << "      Affected variants: " << globalstats.MunphasedRefOnly << "\n";
            yamlinfo << "      Affected variants fraction: " << (globalstats.MunphasedRefOnly / (double) globalstats.Mref) << "\n";
            StatusFile::addWarning("Reference for imputation contains unphased genotypes (set to random phase).<br>\n"
                    "   Variants with unphased data: " + to_string(globalstats.MunphasedRefOnly)
                    + " Fraction: " + to_string(globalstats.MunphasedRefOnly / (double) globalstats.Mref) + "<br>\n"
                    "   Unphased genotypes: " + to_string(globalstats.GunphasedRefOnly)
                    + " Fraction: " + to_string((globalstats.GunphasedRefOnly / (double) globalstats.Mref) / Nref));
        }
    }

    double missingrate = (globalstats.GmissingTarget / (double) globalstats.M) / Ntarget;
    StatusFile::addInfo("Av. missing rate in target genotypes: " + to_string(missingrate) + "</p>", false);
    yamlinfo << "    Missing rate in target genotypes: " << missingrate << "\n";

    if (doImputation) {
        double conflictrate = 0.0;
        if (conflictcnt) {
            size_t num_sites_no_ov = site_offsets.back() - num_sites_per_file[0]; // need to remove the first overlap to get the number of sites without overlaps
            size_t total_haps = num_sites_no_ov * ((getNTarget()-getNHaploidsTgt())*2 + getNHaploidsTgt());
            conflictrate = conflictcnt/(double)total_haps;
            StatusFile::addInfo("<p class='pinfo'>Imputation conflict rate (relation of haplotypes with no overlap and dosage set to allele frequency to total haplotypes): "
                    + to_string(conflictrate) + "</p>", false);
        } else {
            StatusFile::addInfo("<p class='pinfo'>No imputation conflicts.</p>", false);
        }
        yamlinfo << "    Imputation conflict rate: " << conflictrate << "\n";

        if (improveCalls) {
            cout << "\nResults from using --improveCalls:" << endl;
            cout << "  Dosage improvements (no call change): " << improvecntnochange << endl;
            cout << "  Average dosage improvement (no call change): " << (totalimprovenochange / improvecntnochange) << endl;
            cout << "  Dosage improvements (phase switch): " << improvecnthet2het << endl;
            cout << "  Average dosage improvement (phase switch): " << (totalimprovehet2het / improvecnthet2het) << endl;
            cout << "  Dosage improvements (gt switch hom2hom): " << improvecnthom2hom << endl;
            cout << "  Average dosage improvement (gt switch hom2hom): " << (totalimprovehom2hom / improvecnthom2hom) << endl;
            cout << "  Dosage improvements (gt switch hom2het): " << improvecnthom2het << endl;
            cout << "  Average dosage improvement (gt switch hom2het): " << (totalimprovehom2het / improvecnthom2het) << endl;
            cout << "  Dosage improvements (gt switch het2hom): " << improvecnthet2hom << endl;
            cout << "  Average dosage improvement (gt switch het2hom): " << (totalimprovehet2hom / improvecnthet2hom) << endl;
        }
    }

    // write YAML
    if (yaml)
        StatusFile::addInfoYAML("Summary", yamlinfo.str());

    StatusFile::clearContext();
}
