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

#ifndef TARGETIMP_H_
#define TARGETIMP_H_

#include "VCFData.h"
#include "PBWT.h"
#include "Datatypes.h"
#include "RingBuffer.h"

//// DEBUG
//#define COUNTSMMATCHES

using namespace std;

class TargetImp {

public:
    TargetImp(size_t hap_id, const VCFData& vcfdata, const PBWT &pbwt, const BooleanVector &phasedTarget, unsigned num_blocks, int maxerrors, bool imputeCalls, bool debug);

    // calculates the set-maximal matches for the complete target as preparation for imputation
    // and forwards the block states correspondingly (if there is more than one block)
    void calcSMMatches();

    // imputes the next nsites variants according to pre-calculated set-maximal matches and stores the results in imputedTargetBunch and imputedDosageBunch
    void imputeBunch(unsigned block, size_t nsites, BooleanVector &imputedTargetBunch, vector<float> &imputedDosageBunch);

    // DEBUG
#ifdef COUNTSMMATCHES
    size_t num_matches = 0;
    size_t num_sites = 0;
    double av_match_length = 0.0;
    double av_match_height = 0.0;
    size_t reducable_sites = 0;
    void printSMMatches();
#endif
    // __DEBUG

private:

    // class to store the information of a set-maximal match
    class Match {
    public:
        Match() : start(0), end(0), refstart(0), refend(0) {}
        Match(size_t start_, size_t end_, size_t refstart_, size_t refend_) : start(start_), end(end_), refstart(refstart_), refend(refend_) {}
        size_t start;    // start position (SNP), inclusive
        size_t end;      // end position (SNP), exclusive
        size_t refstart; // first ref index in PBWT sort order at end position of match (end-1)
        size_t refend;   // last ref index (inclusive) in PBWT sort order at end position of match (end-1)
#ifdef COUNTSMMATCHES
        size_t reducable_length = 0;
#endif
    };

#ifdef COUNTSMMATCHES
    void printSMMatch(const Match &m);
#endif

    // a tree consisting of these nodes stores PBWT intervals in the history when finding set-maximal-matches
    // @depth depth of the tree until this node, represents the length of the path represented by this node
    // @errors number of errors for the path represented by this node
    class HistoryTreeNode {
    public:
        HistoryTreeNode() : left(NULL), right(NULL), depth(0), errors(0) {}
        HistoryTreeNode(PBWTInterval i, size_t d, int e) : left(NULL), right(NULL), interval(i), depth(d), errors(e) {}
        HistoryTreeNode* left;
        HistoryTreeNode* right;
        PBWTInterval interval;
        size_t depth;
        int errors;
    };

    // processes a subtree of the history tree, i.e. it maps the interval stored in each node:
    // - the to0 part will overwrite the current contents in ltree,
    // - the to1 part will create a new node at rtree (so it has to be NULL on entry)
    // m: current cursor position
    // curr_hap: current haplotype at cursor
    // missing: haplotype is missing
    // treedepth: depth of the tree before start of recursion
    // perase: parent will be erased after processing -> mapping not required, just descend to find longest match (and erase the node afterwards)
    // appends potential matches to "matches"
    // returns the depth of the deepest branch going through either ltree or rtree after mapping
    size_t processHistorySubTree(HistoryTreeNode* &ltree, HistoryTreeNode* &rtree, size_t m, bool curr_hap, bool missing, size_t treedepth, bool perase, vector<Match> &matches);

    // forward imputation state for block by nsites
    void forwardBlock(unsigned block, size_t nsites);

    size_t hap_id; // haplotype ID: /2 -> target ID, %2 -> mat/pat
    const VCFData& vcfdata;
    const PBWT& pbwt;
    const BooleanVector &phasedTarget;

    size_t K; // number of reference haplotypes in PBWT
    size_t M; // number of common sites in reference and targets
    size_t Mref; // number sites in full reference

    vector<Match> sm_matches;
    size_t num_blocks;

    // for the current imputation status (for each block):
    vector<bool> beforefirstcommons; // indicates that we are in the segment before the first common site
    vector<size_t> currms; // last common site (frame index) if beforefirstcommon == false
    vector<size_t> mmaps; // NEXT common site mapped to index in full reference
    vector<size_t> idx0s; // points to the first match probably including mref (i.e. the end > m, should be the case for the first match)
    const RingBuffer<size_t> &misspos;
    vector<size_t> miss_idxs;
    vector<size_t> nextmisss;
    vector<size_t> mrefs; // curr position (in reference) for imputation

    int setmaxerrors; // how many errors are allowed in a set-maximal match?
    bool imputeCalls; // should imputation also be performed explicitly at call sites? since these sites are used as anchors, the result is 100% clear, besides setmaxerrors is > 0!

    bool debug;
};



#endif /* TARGETIMP_H_ */
