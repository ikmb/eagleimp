/*
 *    Copyright (C) 2018-2021 by Lars Wienbrandt and Jan Christian Kässens,
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
#include <cmath>
#include <thread>

#include <boost/program_options.hpp>

#include "Args.h"
#include "VCFData.h"
#include "StatusFile.h"
#include <version.h>

namespace bpo = boost::program_options;

using namespace std;
using namespace bpo;

/**
 * @brief Constructs, parses and verifies all command-line arguments
 *
 * This function constructs the Args object and parses all given command-line
 * arguments. If unknown arguments and/or missing obligatory arguments are
 * detected, this function does not return and instead prints an appropriate
 * help message and calls exit(), executing all exit handlers registered
 * up to the parseArgs call.
 *
 * @param argc argument count
 * @param argv argument vector
 * @return an rvalue reference of a new Args object containing all defined or defaulted CLI arguments
 */
/*static*/ Args Args::parseArgs(int argc, char *argv[]) {
    Args args {argc, argv};

    if (argc <= 1 || args.count("help")) {
        args.printHelp(argv[0], cout);
        exit(EXIT_SUCCESS);
    }

    if (args.count("version")) {
        args.printVersion();
        exit(EXIT_SUCCESS);
    }

    // set variables
    args.parseVars();

    if (!args.count("target") && !args.count("makeQref")) {
        StatusFile::addError("No VCF target genotype file specified. Try --help.");
        exit(EXIT_FAILURE);
    }

    if (!args.count("ref")) {
        StatusFile::addError("No VCF or Qref reference file specified. Try --help.");
        exit(EXIT_FAILURE);
    }
    size_t qpos = args.ref.rfind(".qref");
    args.isQRef = qpos != string::npos && qpos == args.ref.size()-5;

    if (args.count("makeQref") && args.isQRef) {
        StatusFile::addError("Qref creation requires a VCF reference. Try --help.");
        exit(EXIT_FAILURE);
    }

    // checked later due to special case with chrY where no phasing and thus no genetic map is required
//    if (!args.count("geneticMap") && !args.count("makeQref") && !args.count("skipPhasing")) {
//        StatusFile::addError("No genetic map file specified. Try --help.");
//        exit(EXIT_FAILURE);
//    }

    if (!args.count("output")) { // take the base of the target file as output prefix
        if (args.vcfTarget.size() >= 7 && args.vcfTarget.substr(args.vcfTarget.size()-7).compare(".vcf.gz") == 0) // special case for .vcf.gz
            args.outPrefix = args.vcfTarget.substr(0, args.vcfTarget.size() - 7);
        else // standard case: simply leave out file suffix starting with '.'
            args.outPrefix = args.vcfTarget.substr(0, args.vcfTarget.rfind('.'));
    }

    // set write mode according to desired output mode
    if (!args.vcfOutFormat.compare("b"))
        args.writeMode = string("wb");
    else if (!args.vcfOutFormat.compare("u"))
        args.writeMode = string("wbu");
    else if (!args.vcfOutFormat.compare("z"))
        args.writeMode = string("wz");
    else if (!args.vcfOutFormat.compare("v"))
        args.writeMode = string("w");
    else {
        StatusFile::addError("Invalid VCF output format. Please specify either of b|u|z|v.");
        exit(EXIT_FAILURE);
    }

    // compare filename arguments to be unique
    args.checkFilenames();

    if (args.num_threads > std::thread::hardware_concurrency()) {
        StatusFile::addWarning("The number of threads is higher than the number of available system threads. Reduced number of threads to " + to_string(std::thread::hardware_concurrency()) + ".");
        args.num_threads = std::thread::hardware_concurrency();
    }

    if (args.num_files > args.num_threads) {
        StatusFile::addError("The number of temporary files must not exceed the number of worker threads.");
        exit(EXIT_FAILURE);
    }

    if (args.bpEnd >= MAX_CSI_COOR) {
        StatusFile::addError("End position must not exceed maximum indexable position " + to_string(MAX_CSI_COOR));
        exit(EXIT_FAILURE);
    }

    if (args.bpEnd < args.bpStart) {
        StatusFile::addError("Invalid region! bpEnd < bpStart");
        exit(EXIT_FAILURE);
    }

    if (args.K < 2) {
        StatusFile::addError("K must at least be 2, though we do not recommend reducing K below the default setting.");
        exit(EXIT_FAILURE);
    }

    return args;
}

ostream &operator<<(ostream &out, const Args &args) {
    variables_map::const_iterator it;

    long name_width = 0;
    long value_width = 0;
    long orig_width = out.width();

    // collect field widths
    for(it = args.vars.begin(); it != args.vars.end(); ++it) {
        long this_width = static_cast<long>(it->first.length());

        if(this_width > name_width) {
            name_width = this_width;
        }

        this_width = 0;

        if(it->second.value().type() == typeid(string)) {
            this_width = static_cast<long>(it->second.as<string>().length());
        }

        if(it->second.value().type() == typeid(int)) {
            this_width = static_cast<long>(log10(it->second.as<int>()));
        }

        if(this_width > value_width) {
            value_width = this_width;
        }
    }

    // dump option values
    out.setf(ios::adjustfield, ios::left);

    for(it = args.vars.begin(); it != args.vars.end(); ++it) {
        out.width(name_width + 2); // more space for the colon
        out << (it->first + ":") << endl;

        out.width(value_width);

        if(it->second.value().type() == typeid(string)) {
            out << it->second.as<string>() << endl;
        } else if(it->second.value().type() == typeid(int)) {
            out << it->second.as<int>() << endl;
        } else {
            out << "(unknown)" << endl;
        }
    }

    resetiosflags(ios::adjustfield);
    out.width(orig_width);

    return out;
}

Args::Args(int argc, char *argv[]) :
    opts_regular("Program options"),
	opts_regular_hybrid("Additional program options for hybrid architecture"),
    opts_region("Region selection options, the chromosome is always set automatically from the filename of the reference file"),
	opts_phasing("Phasing algorithm options"),
	opts_imputation("Imputation algorithm options"),
	opts_hidden("Hidden options (only visible in debug mode)"),
	opts_hidden_hybrid("Additional hidden options for hybrid architecture")
	{

    opts_regular.add_options()
    ("help,h", "produces this help message and exits")
    ("version,v", "prints version information and exits")
    ("threads,t", value<unsigned>(&num_threads)->default_value(std::thread::hardware_concurrency()), "number of threads to use for phasing and imputation")
    ("output,o", value<string>(&outPrefix), "output file prefix (optional)")
    ("geneticMap", value<string>(&genMapFile), "HapMap genetic map")
    ("ref", value<string>(&ref), "tabix-indexed compressed VCF/BCF file or Qref file for reference haplotypes")
    ("target", value<string>(&vcfTarget), "tabix-indexed compressed VCF/BCF file for target genotypes")
    ("vcfExclude", value<string>(&vcfExclude), "tabix-indexed compressed VCF/BCF file containing variants to exclude from phasing")
    ("vcfOutFormat", value<string>(&vcfOutFormat)->default_value("z"), "b|u|z|v: compressed BCF (b), uncomp BCF (u), compressed VCF (z), uncomp VCF (v)")
    ("excludeMultiAllRef", "exclude multi-allelic reference markers for phasing as well as for imputation")
    ("allowRefAltSwap", "allow swapping of REF/ALT in target vs. reference VCF (the output will match reference alleles)")
    ("allowStrandFlip", "allow strand flips, i.e. A/T and C/G swaps in target vs. reference (precedence for allowRefAltSwap on A/T and C/G variants!) (the output will match reference alleles)")
    ("skipPhasing", "skip phasing step prior to imputation, instead the input target is assumed to be phased")
    ("skipImputation", "skip reference imputation")
    ("makeQref", "create quick reference by dumping the decoded VCF reference, --ref is required with a VCF reference file, all other args are ignored")
    ("maxChunkMemory", value<uint64_t>(&maxChunkMemory)->default_value(16ull), "in GiB, targets that require a pre-estimated amount of memory larger than this will be divided into chunks")
    ("chunkFlankSize", value<uint64_t>(&chunkFlankSize)->default_value(64), "number of target variants that will be used as flanking region to either side of a chunk; this will not extend a previously set region!")
    ("stat", value<string>(&statfile), "file for status output")
    ("lockfile", value<string>(&lockfile), "optional lockfile for exclusive CPU usage (used by CPU phasing and imputation). if not present, no CPU lock will be used.")
    ;

#ifdef USE_AD_FPGA
    opts_regular_hybrid.add_options()
    ("fpgas", value<vector<unsigned>>(&fpgas)->multitoken()->default_value(vector<unsigned>(),""), "activate FPGA acceleration and restrict to these FPGAs (separate indices by whitespace)")
    ;
#endif
#ifdef USE_CUDA_GPU
    opts_regular_hybrid.add_options()
    ("gpus", value<vector<unsigned>>(&gpus)->multitoken()->default_value(vector<unsigned>(),""), "activate GPU acceleration and restrict to these GPUs (separate indices by whitespace)")
    ;
#endif

    opts_region.add_options()
    // the chromosome is always set automatically
    //("chrom", value<string>(&chromStr)->default_value("0"), "chromosome to analyze (if input has many, 0=first (either read from target or reference!))")
    ("bpStart", value<int64_t>(&bpStart)->default_value(0), "minimum base pair position to analyze (1-based)")
    ("bpEnd", value<int64_t>(&bpEnd)->default_value(MAX_CSI_COOR-1), "maximum base pair position to analyze (inclusive, 1-based)")
    ("bpFlanking", value<int64_t>(&bpFlanking)->default_value(0), "flanking region to use during phasing/imputation but discard in output")
    ;

    opts_phasing.add_options()
    ("doPrePhasing", "explicitly enable prior pre-phasing step. Pre-phasing is disabled per default as benchmarks showed no benefits, neither in performance nor in quality.")
    ("noRevPhasing", "disable posterior reverse phasing step (not recommended!)")
    ("impMissing", "enable imputation of missing target genotypes (. or ./.) during phasing (missings are filled according to the ref panel AF), otherwise missings are imputed during imputation after phasing")
    ("outputUnphased", "output unphased sites in phased output file (including excluded variants), anyway unphased sites are excluded in imputation output!")
    ("Kpbwt,K", value<size_t>(&K)->default_value(10000), "number of conditioning haplotypes (0=maximum - use with caution!)")
    ("phaseIters,i", value<uint32_t>(&iters)->default_value(0), "number of PBWT phasing iterations (0=auto)")
    ("expectIBDcM", value<fp_type>(&expectIBDcM)->default_value(2.0), "expected length of haplotype copying (cM)")
    ("histFactor", value<fp_type>(&hist)->default_value(0.0), "history length multiplier (0=auto)")
    ("genoErrProb", value<fp_type>(&pErr)->default_value(0.003), "estimated genotype error probability")
    ;

    opts_imputation.add_options()
    ("imputeInfo", value<string>(&impInfoStr)->default_value("a"), "string indicating desired imputation information for each genotype: 'a' allele dosage, 'g' genotype dosage, 'p' genotype probabilities; argument may contain any combination or none")
    ("imputeFilter", value<double>(&impFilter)->default_value(0.0), "r2 score threshold for imputation output")
    ("outputPhasedFile", "output phased sites in extra file if reference imputation is enabled, otherwise phased sites are written anyway")
    ("noMissingIDs", "convert missing variant IDs to chr:pos:ref:alt in imputation output or Qref, this option is automatically enabled when creating a Qref")
    ;

    opts_hidden.add_options()
    ("cMmaxsplit", value<float>(&cMmaxsplit)->default_value(0.125), "max run of homozygous genotypes for phasing in cM")
//    ("minibdA", value<size_t>(&min_ibd_a)->default_value(40), "minimum number of consistent sites for IBD A condition (zero disables IBD check)")
//    ("minibdB", value<size_t>(&min_ibd_b)->default_value(20), "minimum number of consistent sites for IBD B condition (must be <= minibdA - zero disables IBD check)")
    ("hsize", value<size_t>(&hsize)->default_value(100), "history size, must not be larger than 128 but larger than delta+1")
    ("hsize_fast", value<size_t>(&hsize_fast)->default_value(30), "history size in fast pre-phasing, must not be larger than 128, but larger than delta_fast+2")
    ("beamwidth", value<size_t>(&beamwidth)->default_value(128), "beam width")
    ("beamwidth_fast", value<size_t>(&beamwidth_fast)->default_value(64), "beam width in fast pre-phasing")
    ("delta", value<size_t>(&delta)->default_value(24), "cursor is advanced delta sites before actual site is processed")
    ("delta_fast", value<size_t>(&delta_fast)->default_value(10), "delta in fast pre-phasing")
    ("minfixsize", value<float>(&minfixsize)->default_value(0.5), "relative amount of sites that at least have to be constrained in fast pre-phasing")
    ("maxfixsize", value<float>(&maxfixsize)->default_value(0.9), "relative amount of sites that should be constrained at most in fast pre-phasing (if not >=maxfixthr confident)")
    ("minfixthr", value<float>(&minfixthr)->default_value(0.99), "constrained in pre-phasing if confidence is greater or equal in the area between minfixsize and maxfixsize")
    ("maxfixthr", value<float>(&maxfixthr)->default_value(1.0), "constrained in pre-phasing if confidence is greater or equal, no matter what!")
    ("bplimitSq", "use pErr^2 as maximum probability deviation in beam instead of pErr")
    ("imputeBunchSize", value<size_t>(&ibunchsize)->default_value(10000), "number of fields (x 1000) to be imputed at a time (variants x samples x 1000)")
    ("setmaxerrors", value<int>(&setmaxerrors)->default_value(0), " number of maximal tolerated errors in a set-maximal match")
    ("overwriteCalls", "imputation is allowed to change genotype calls from the input, it will do so whenever the imputation result is different to the original call")
    ("improveCalls", "as in --overwriteCalls imputation is allowed to change genotype calls, but only if the dosage will be improved")
    ("numTmpFiles", value<unsigned>(&num_files)->default_value(0), "number of temporary files for imputation output (0=auto)")
    ("overrideChunks", value<int>(&overrideChunks)->default_value(0), "disable automatic chunk division and set to provided number (0=auto enabled)")
    ("skipHeader", "skip writing of VCF header for output files")
    ("deterministic", "disable random generator such that the output is deterministic for each run")
    ("debug", "produce lots of debug output")
    ;

#ifdef USE_CUDA_GPU
    opts_hidden_hybrid.add_options()
    ("buffer-size-GPU", value<size_t>()->default_value(2ull*1024*1024*1024), "Size for transmission buffers (FPGA->host) in bytes.")
    ("buffers-GPU", value<unsigned>()->default_value(0), "Number of transmission buffers (GPU->host) to keep around. (0=automatic)")
    ;
#endif
#ifdef USE_AD_FPGA
    opts_hidden_hybrid.add_options()
    ("buffer-size-FPGA", value<size_t>()->default_value(128*1024*1024), "Size for transmission buffers (FPGA->host) in bytes.")
    ("buffers-FPGA", value<unsigned>()->default_value(0), "Number of transmission buffers (FPGA->host) to keep around. (0=automatic)")
    ("timeout", value<unsigned long>(&timeout)->default_value(600000), "Timeout for FPGA transmissions (in ms)")
    ;
#endif

    parse(argc, argv);
}

void Args::parse(int argc, char *argv[]) {
    bpo::options_description all_options;

    // combine all options
    all_options.add(opts_regular);
#if defined USE_CUDA_GPU || defined USE_AD_FPGA
    all_options.add(opts_regular_hybrid);
#endif
    all_options.add(opts_region).add(opts_phasing).add(opts_imputation).add(opts_hidden);
#if defined USE_CUDA_GPU || defined USE_AD_FPGA
    all_options.add(opts_hidden_hybrid);
#endif

    // do the actual parsing
    store(command_line_parser(argc, argv).options(all_options).run(), vars);
    notify(vars);

}

void Args::parseVars() {

    if (!statfile.empty()) {
        bool ok = true;
        size_t pos = statfile.rfind(".vcf.gz");
        if (pos != string::npos && pos == statfile.size()-7)
            ok = false;
        pos = statfile.rfind(".vcf");
        if (pos != string::npos && pos == statfile.size()-4)
            ok = false;
        pos = statfile.rfind(".bcf");
        if (pos != string::npos && pos == statfile.size()-4)
            ok = false;
        pos = statfile.rfind(".qref");
        if (pos != string::npos && pos == statfile.size()-5)
            ok = false;
        if(!ok) {
            // Due to a potential overwrite of significant data, this error message is not written in the status file, but printed only.
            cerr << "ERROR: Status file must not end with .vcf.gz/.vcf/.bcf/.qref" << endl;
            exit(EXIT_FAILURE);
        }
        StatusFile::updateFile(statfile);
    }

    // set bools
    if (vars.count("impMissing"))
        impMissing = true;
    if (vars.count("allowRefAltSwap"))
        allowRefAltSwap = true;
    if (vars.count("allowStrandFlip"))
        allowStrandFlip = true;
    if (vars.count("excludeMultiAllRef"))
        excludeMultiAllRef = true;
    if (vars.count("outputUnphased"))
        outputUnphased = true;
    if (vars.count("noMissingIDs"))
        noMissingIDs = true;
    if (vars.count("outputPhasedFile"))
        outputPhased = true;
    if (vars.count("skipPhasing"))
        skipPhasing = true;
    if (vars.count("skipImputation"))
        skipImputation = true;
    if (vars.count("doPrePhasing"))
        doPrePhasing = true;
    if (vars.count("noRevPhasing"))
        noRevPhasing = true;
    if (vars.count("bplimitSq"))
        bplimitSq = true;
    if (vars.count("overwriteCalls"))
        overwriteCalls = true;
    if (vars.count("improveCalls"))
        improveCalls = true;

    // parse imputation writeout information
    if (impInfoStr.find('a') != string::npos)
        writeADosage = true;
    if (impInfoStr.find('g') != string::npos)
        writeGDosage = true;
    if (impInfoStr.find('p') != string::npos)
        writeProb = true;

    // automatically define number of temporary files for imputation output
    if (num_files == 0)
        num_files = max(1u, num_threads/4); // empirically

    // dump reference data
    if (vars.count("makeQref")) {
        createQuickRef = true;
        // there will be no phasing or imputation!
        skipPhasing = true;
        skipImputation = true;
        // missing IDs are automatically converted
        noMissingIDs = true;
    }

    // set K to max, if desired
    if (K == 0)
        K = ~0ull;

    // skip VCF header
    if (vars.count("skipHeader"))
        skipHeader = true;

    // deterministic
    if (vars.count("deterministic"))
        deterministic = true;

    // debug
    if (vars.count("debug"))
        debug = true;

}

// compare filename arguments to be unique
inline void Args::checkFilenames() {

    // check input files
    {
        if (!vcfTarget.compare(ref)) {
            StatusFile::addError("Target and reference cannot be the same file.");
            exit(EXIT_FAILURE);
        }
        if (!vcfExclude.empty() && (!vcfExclude.compare(ref) || !vcfExclude.compare(vcfTarget))) {
            StatusFile::addError("Exclude file cannot be equal to target or reference.");
            exit(EXIT_FAILURE);
        }
        // other crazy combinations with equal files will fail anyway due to wrong formats (e.g. if a VCF file equals the genetic map file)
    }

    // check status file
    if (!statfile.empty()) {
        if(!statfile.compare(outPrefix)) {
            StatusFile::addError("Status file and output prefix cannot be equal.");
            exit(EXIT_FAILURE);
        }
        if(!statfile.compare(lockfile)) {
            StatusFile::addError("Status file and lock file cannot be equal.");
            exit(EXIT_FAILURE);
        }
        // status file was already checked to not have a file suffix as common input files
    }

    // check lock file
    if (!lockfile.empty()) {
        {
            bool ok = true;
            size_t pos = lockfile.rfind(".vcf.gz");
            if (pos == lockfile.size()-7)
                ok = false;
            pos = lockfile.rfind(".vcf");
            if (pos == lockfile.size()-4)
                ok = false;
            pos = lockfile.rfind(".bcf");
            if (pos == lockfile.size()-4)
                ok = false;
            pos = lockfile.rfind(".qref");
            if (pos == lockfile.size()-5)
                ok = false;
            pos = lockfile.rfind(".varinfo");
            if (pos == lockfile.size()-8)
                ok = false;
            pos = lockfile.rfind(".info");
            if (pos == lockfile.size()-5)
                ok = false;
            pos = lockfile.rfind(".warning");
            if (pos == lockfile.size()-8)
                ok = false;
            pos = lockfile.rfind(".error");
            if (pos == lockfile.size()-6)
                ok = false;
            if(!ok) {
                StatusFile::addError("Lock file must not end with .vcf.gz/.vcf/.bcf/.qref/.varinfo/.info/.warning/.error");
                exit(EXIT_FAILURE);
            }
        }
    }
}

bool Args::isDefined(const string &optname) const {
    bool found = false;
    found = !!this->opts_regular.find_nothrow(optname, false); // return null (-> false) if option has not been found
    found |= !!this->opts_hidden.find_nothrow(optname, false);
    found |= !!this->opts_region.find_nothrow(optname, false);
    found |= !!this->opts_phasing.find_nothrow(optname, false);
    found |= !!this->opts_imputation.find_nothrow(optname, false);
#if defined USE_CUDA_GPU || defined USE_AD_FPGA
    	found |= !!this->opts_regular_hybrid.find_nothrow(optname, false);
    	found |= !!this->opts_hidden_hybrid.find_nothrow(optname, false);
#endif
    return found;
}

void Args::printHelp(const string &progname, ostream &out) const {
    out << "Usage: " << progname << " [options]" << endl << endl;
    out << opts_regular << endl;
#if defined USE_CUDA_GPU || defined USE_AD_FPGA
    out << opts_regular_hybrid << endl;
#endif
    out << opts_region << endl;
    out << opts_phasing << endl;
    out << opts_imputation << endl;
    if (debug) {
        out << opts_hidden << endl;
#if defined USE_CUDA_GPU || defined USE_AD_FPGA
        out << opts_hidden_hybrid << endl;
#endif
    }
    out << endl;

    printVersion();
}

/* static */
void Args::printVersion() {
    cout << "This is version " << prog_version << ", compiled on " << prog_timestamp << endl;
    cout << "Send bugs to " << prog_bugaddress << endl;
}
