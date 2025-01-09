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

#ifndef ARGS_H_
#define ARGS_H_

#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "Datatypes.h"

namespace bpo = boost::program_options;

using namespace std;

/**
 * Class for storing and retrieving command-line arguments.
 */
class Args {
public:

	static Args parseArgs(int argc, char *argv[]);

    /**
     * Returns the argument with the given (long) name. The template argument
     * specifies the return type you want the parameter casted to.
     * @param name argument name
     * @return the variable value
     */
    template<typename T> T get(const std::string &name) const {
        auto where = vars.find(name);
        if(where == std::end(vars)) {
            if(!isDefined(name))
                throw std::invalid_argument("Option undefined: " + name + " (This is a bug)");
            else
                throw std::out_of_range("Option has not been specified and does not have a default value associated: " + name);
        }
        return where->second.as<T>();
    }

    /**
     * Counts the number of argument occurences. Mainly useful for boolean switches
     * and/or counting flags, such as verbosity or debug mode.
     * @param name argument name
     * @return argument value
     */
    unsigned int count(const std::string &name) const {
        if(!isDefined(name))
            throw std::invalid_argument("Option undefined: " + name + " (This is a bug)");
        return vars.count(name);
    }

    bool operator()(const std::string &name) const {
        return count(name) > 0;
    }

    /**
     * Prints a help message.
     * @param progname the program name, usually argv[0]
     * @param out output stream, e.g. cout
     */
    void printHelp(const std::string &progname, std::ostream &out) const;
    static void printVersion();

    Args(Args&& other) = default;

    unsigned num_threads;
    string outPrefix;
    string genMapFile;
    string ref;
    bool isQRef;
    string vcfTarget;
    string vcfExclude;
    string writeMode;
    bool impMissing = false;
    bool allowRefAltSwap = false;
    bool allowStrandFlip = false;
    bool excludeMultiAllRef = false;
    bool outputUnphased = false;

    bool skipImputation = false;
    float impR2filter;
    float impMAFfilter;
    bool noMissingIDs = false;

    bool outputPhased = false;
    bool skipPhasing = false;

    bool writeADosage = false;
    bool writeGDosage = false;
    bool writeProb = false;

    vector<unsigned> fpgas;
    vector<unsigned> gpus;

//    int chrom;
    // htslib uses int64_t...
    int64_t bpStart;
    int64_t bpEnd;
    int64_t bpFlanking;
    uint64_t maxChunkMemory;
    uint64_t chunkFlankSize;

    size_t K;
    uint32_t iters;
    fp_type expectIBDcM;
    fp_type hist;
    fp_type pErr;

    float  cMmaxsplit;
//    size_t min_ibd_a;
//    size_t min_ibd_b;
    bool   doPrePhasing = false;
    bool   noRevPhasing = false;
    size_t hsize;
    size_t hsize_fast;
    size_t beamwidth;
    size_t beamwidth_fast;
    size_t delta;
    size_t delta_fast;
    float  minfixsize;
    float  maxfixsize;
    float  minfixthr;
    float  maxfixthr;
    bool   bplimitSq = false;
    int setmaxerrors;
    bool overwriteCalls = false;
    bool improveCalls = false;
    size_t ibunchsize;
    unsigned num_files;

    string statfile;
    string lockfile;

    bool createQuickRef = false;
    bool skipHeader = false;

    int overrideChunks;
    bool deterministic = false;
    bool debug = false;
    bool yaml = false;

    unsigned long timeout;

protected:
    /** Constructs the arguments list and adds all defined options */
    Args();
    Args(Args const &);
    void operator=(Args const &);

    void parse(int argc, char *argv[]);
    void parseVars();
    bool isDefined(const std::string &optname) const;

    bpo::options_description opts_regular;        /**< regular options, shown on help output */
    bpo::options_description opts_regular_hybrid; /**< regular options, shown on help output */
    bpo::options_description opts_region;         /**< regular region select options, shown on help output */
    bpo::options_description opts_phasing;        /**< regular phasing algorithm options, shown on help output */
    bpo::options_description opts_imputation;     /**< regular imputation algorithm options, shown on help output */
    bpo::options_description opts_hidden;         /**< hidden options */
    bpo::options_description opts_hidden_hybrid;  /**< hidden options */

    bpo::variables_map vars;    /**< this is where the option values go */

    /** dump all options to the given stream */
    friend std::ostream &operator<<(std::ostream &out, const Args &args);

private:
    /** parses the main() options into the variable map */
    Args(int argc, char *argv[]);

    inline void checkFilenames();

    string vcfOutFormat;
    string impInfoStr;
};

#endif /* ARGS_H_ */
