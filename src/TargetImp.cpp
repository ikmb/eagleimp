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
#include "PBWTInterval.h"

#include "TargetImp.h"

static const float DOSAGE_CONFLICT = -1.0;

TargetImp::TargetImp(size_t hap_id_, const VCFData& vcfdata_, const PBWT &pbwt_, const BooleanVector &phasedTarget_, unsigned num_blocks_, int setmaxerrors_, bool imputeCalls_, bool debug_) :
    hap_id(hap_id_),
    vcfdata(vcfdata_),
    pbwt(pbwt_),
    phasedTarget(phasedTarget_),
    K(pbwt_.getK()),
    M(vcfdata.getNSNPs()),
    Mref(vcfdata.getNSNPsFullRef()),
    num_blocks(num_blocks_),
    misspos(vcfdata_.getTgtMissPos()[vcfdata_.getHaploidsTgtMap()[hap_id_]/2]),
    setmaxerrors(setmaxerrors_),
    imputeCalls(imputeCalls_),
    debug(debug_) {
    sm_matches.reserve(M/8); // empirically, works with any value

    // initialize the imputation status:
    beforefirstcommons.clear();
    currms.clear();
    mmaps.clear();
    idx0s.clear();
    miss_idxs.clear();
    nextmisss.clear();
    mrefs.clear();

    beforefirstcommons.resize(num_blocks, true); // indicates that we are in the segment before the first common site
    currms.resize(num_blocks, 0); // last common site (frame index) if beforefirstcommon == false
    mmaps.resize(num_blocks, vcfdata.getSNPIndexToFullRef(0)); // NEXT common site mapped to index in full reference
    idx0s.resize(num_blocks, 0); // points to the first match probably including mref (i.e. the end > m, should be the case for the first match)

    miss_idxs.resize(num_blocks, 0);
    nextmisss.resize(num_blocks, ~0ull);
    for (unsigned block = 0; block < num_blocks; block++) {
        if (miss_idxs[block] < misspos.size()) {
            nextmisss[block] = misspos[miss_idxs[block]];
            miss_idxs[block]++;
        }
    }
    mrefs.resize(num_blocks, 0);

    // after calculating the sm matches the states have to be forwarded to the corresponding beginning of their block with prepareBlocks()!!!
}

// find all set maximal matches for this target
// the set-maximal matches will be sorted by their "end" field
// AND according to the definition of a set-maximal match, as a result, the matches will also be sorted by their "start" field
// (if not, that would indicate that one match includes another completely, which is not possible by definition)
void TargetImp::calcSMMatches() {

    // We check the missingness rate of the target here.
    // If the missingness rate is above 50%, we do not report any set-maximal matches.
    // As this is usually called in parallel for several targets, we should not generate a warning here.
    if (misspos.size() > M/2)
        return;

    // keep a history of the last intervals in the reference PBWT matching the target
    // from the current position until beginning of first empty interval.
    // the container contains information on the interval,
    // the length of the subsequence represented by the interval
    // and the number of tolerated errors in the entry
    HistoryTreeNode *history_root = new HistoryTreeNode(PBWTInterval(0,K-1),0,0);

    // get position of first missing hap
    const auto &curr_misspos = misspos;
    size_t missidx = 0;
    size_t nextmiss = ~0ull;
    if (missidx < curr_misspos.size()) {
        nextmiss = misspos[missidx];
        missidx++;
    }

    // iterate over all target sites
    size_t treedepth = 0;
    for (size_t m = 0; m < M; m++) {
        bool curr_hap = phasedTarget[m];
        bool missing = m == nextmiss;
        if (missing) {
            // find next missing
            if (missidx >= curr_misspos.size())
                nextmiss = ~0ull;
            else {
                nextmiss = curr_misspos[missidx];
                missidx++;
            }
        }

        // process the history tree:
        // insert complete interval as new tree root, size is 0
        HistoryTreeNode *newroot = new HistoryTreeNode(PBWTInterval(0,K-1),0,0);
        newroot->left = history_root;
        history_root = newroot;
        vector<Match> newmatches;
        treedepth = processHistorySubTree(history_root->left, history_root->right, m, curr_hap, missing, treedepth, false, newmatches);

//        // DEBUG
//        if (newmatches.size() > 0)
//            cerr << "  New matches: m: " << m << " cnt0: " << pbwt.getCount0(m) << " hap: " << (missing ? "?" : (curr_hap ? "1" : "0")) << " treedepth: " << treedepth << endl;
//        for (const auto &match : newmatches) {
//           cerr << "  [" << match.start << "," << match.end << ")(" << (match.end-match.start) << ") x [" << match.refstart << "," << match.refend << "](" << (match.refend - match.refstart + 1) << ")" << endl;
//        }
//        // __DEBUG

        // The reported matches all have the same size, since they are reported only from the maximum depth level.
        // Anyway, the matches will only be set-maximal if there is no longer path still alive in the tree.
        if (newmatches.size() && newmatches[0].end - newmatches[0].start >= treedepth) {
            for (auto &nm : newmatches) {
                sm_matches.push_back(move(nm));
            }
        }

    } // END for m

    // From the remaining history tree, the last items form the matches including the last site.
    // We get these matches by processing the tree with the erase flag set. This will not do any tree modifications but
    // erasing all nodes and report the matches from the leafs.
    vector<Match> newmatches;
    treedepth = processHistorySubTree(history_root, history_root, M, false, false, treedepth, true, newmatches); // NOTE: the right tree won't be touched. the compiler just needs a valid pointer reference here

//    // DEBUG
//    if (newmatches.size() > 0)
//        cerr << "  Last new matches: m: " << M << " treedepth: " << treedepth << endl;
//    for (const auto &match : newmatches) {
//       cerr << "  [" << match.start << "," << match.end << ")(" << (match.end-match.start) << ") x [" << match.refstart << "," << match.refend << "](" << (match.refend - match.refstart + 1) << ")" << endl;
//    }
//    // __DEBUG

    // The reported matches all have the same size, since they are reported only from the maximum depth level.
    // Since the tree is deleted now, all matches are set-maximal.

    for (auto &nm : newmatches) {
        sm_matches.push_back(move(nm));
    }

    // DEBUG
#ifdef COUNTSMMATCHES
    {
        num_sites = 0;
        size_t mlsum = 0, mhsum = 0;
        for (Match &match : sm_matches) {
            size_t match_length = match.end - match.start;
            size_t match_height = match.refend - match.refstart + 1;
            mlsum += match_length;
            mhsum += match_height;

            // check for which sites quick imputation is possible (if AF is very low)
            size_t mref_start = vcfdata.getSNPIndexToFullRef(match.start);
            size_t mref_end = match.end == vcfdata.getNSNPs() ? vcfdata.getNSNPsFullRef() : vcfdata.getSNPIndexToFullRef(match.end);
            size_t mref_length = mref_end - mref_start;
            size_t m_height = match.refend - match.refstart + 1;
            num_sites += mref_length * m_height;
            size_t nref = vcfdata.getNReferenceHaps();
            size_t mref_reducable_length = 0; // reducable length (sites that can be quickly imputed)
            for (size_t mr = mref_start; mr < mref_end; mr++) {
                float af = vcfdata.getAlleleFreqFullRefRegion(mr);
                float maf = af > 0.5 ? 1-af : af;
                if (maf*nref < 0.1*m_height) {
                    mref_reducable_length++;
                }
            }
            reducable_sites += mref_reducable_length * m_height;
            match.reducable_length = mref_reducable_length;
        }
        num_matches = sm_matches.size();
        av_match_length = mlsum / (double)num_matches;
        av_match_height = mhsum / (double)num_matches;
    }
#endif
    // __DEBUG

//    // DEBUG
//    cerr << "\nSet maximal matches: [start, end) refstart : refend  (#refs)" << endl;
//    int index = 0;
//    for (const auto &match : sm_matches) {
//        cerr << "match " << index << ": [" << match.start << ", " << match.end << ") " << match.refstart << " : " << match.refend << " (" << (match.refend-match.refstart+1) << ")" << endl;
//
//        const size_t context = 3;
//        cerr << "tgt:\t\t";
//        for (size_t m = match.start < context ? 0 : match.start - context; m < (match.end + context > phasedTarget.size() ? phasedTarget.size() : match.end + context); m++) {
//            if (m == match.start || m == match.end)
//                cerr << " | ";
//            bool hap = phasedTarget[m];
//            bool missing = vcfdata.getTargets()[hap_id/2][m] == Genotype::Miss;
//            cerr << (missing ? "?" : (hap ? "1" : "0"));
//        }
//        cerr << endl;
//
//        for (size_t ref = match.refstart; ref <= match.refend; ref++) {
//            size_t absref = pbwt.getSortOrder(match.end)[ref]; // find the corresponding absolute index in the reference
//            cerr << "ref " << absref << ":\t";
//            if (absref < 100)
//                cerr << "\t";
//            int mismatches = 0;
//            int contextmm = 0;
//            for (size_t m = match.start < context ? 0 : match.start - context; m < (match.end + context > phasedTarget.size() ? phasedTarget.size() : match.end + context); m++) {
//                if (m == match.start || m == match.end)
//                    cerr << " | ";
//                size_t mref = vcfdata.getSNPIndexToFullRef(m);
//                bool hap = vcfdata.getFullReferenceT()[mref][absref];
//                bool tgtmissing = vcfdata.getTargets()[hap_id/2][m] == Genotype::Miss;
//                cerr << (hap ? 1 : 0);
//                if (!tgtmissing && m >= match.start && m < match.end && phasedTarget[m] != hap) {
//                    mismatches++;
//                }
//                if (!tgtmissing && (m == match.start-1 || m == match.end) && phasedTarget[m] == hap) { // check only the first left and first right context positions, they must be unequal!
//                    contextmm++;
//                }
//            }
//            if (mismatches > setmaxerrors)
//                cerr << " " << mismatches << " MISMATCHES!";
//            else
//                cerr << " ok (" << mismatches << ")";
//            if (contextmm)
//                cerr << " CONTEXT!";
//            cerr << endl;
//        }
//
//        index++;
//    }
//    // __DEBUG
}

// for the distribution in several blocks (for parallel processing) prepare the start states for each block
// NOTE: the first and the last field in num_sites_per_block should contain the number of vars in the overlap which should be ignored
void TargetImp::prepareBlocks(const vector<size_t>& num_sites_per_block) {

    // forward first block by the number of sites in the overlap
    forwardBlock(0, num_sites_per_block[0]);

    // forward imputation block states
    for (unsigned block = 1; block < num_blocks; block++) {
        // copy state from previous block
        beforefirstcommons[block] = beforefirstcommons[block-1];
        currms[block] = currms[block-1];
        mmaps[block] = mmaps[block-1];
        idx0s[block] = idx0s[block-1];
        miss_idxs[block] = miss_idxs[block-1];
        nextmisss[block] = nextmisss[block-1];
        mrefs[block] = mrefs[block-1];

        // forward state by number of sites in PREVIOUS block (but keeping the extra overlap field in mind!)
        forwardBlock(block, num_sites_per_block[block]);
    }

//    // DEBUG
//    if (hap_id == 0) {
//        cerr << "Mref:\t" << Mref << endl;
//        cerr << "M:\t" << M << endl;
//        cerr << "numsites:";
//        for (size_t n = 0; n < num_blocks+2; n++)
//            cerr << "\t" << num_sites_per_block[n];
//        cerr << endl;
//        cerr << "mrefs:\t\t";
//        for (size_t n = 0; n < num_blocks; n++)
//            cerr << "\t" << mrefs[n];
//        cerr << endl;
//        cerr << "currms:\t\t";
//        for (size_t n = 0; n < num_blocks; n++)
//            cerr << "\t" << currms[n];
//        cerr << endl;
//        cerr << "mmaps:\t\t";
//        for (size_t n = 0; n < num_blocks; n++)
//            cerr << "\t" << mmaps[n];
//        cerr << endl;
//        cerr << "idx0s:\t\t";
//        for (size_t n = 0; n < num_blocks; n++)
//            cerr << "\t" << idx0s[n];
//        cerr << endl;
//        cerr << "missidx:\t";
//        for (size_t n = 0; n < num_blocks; n++)
//            cerr << "\t" << miss_idxs[n];
//        cerr << endl;
//        cerr << "nextmiss:\t";
//        for (size_t n = 0; n < num_blocks; n++)
//            cerr << "\t" << nextmisss[n];
//        cerr << endl;
//    }
//    // __DEBUG
}

#ifdef COUNTSMMATCHES
void TargetImp::printSMMatches() {
    for (const auto &m : sm_matches) {
        printSMMatch(m);

        size_t mref_start = vcfdata.getSNPIndexToFullRef(m.start);
        size_t mref_end = m.end == vcfdata.getNSNPs() ? vcfdata.getNSNPsFullRef() : vcfdata.getSNPIndexToFullRef(m.end);
        size_t mref_length = mref_end - mref_start;
        size_t m_height = m.refend - m.refstart + 1;
        size_t nref = vcfdata.getNReferenceHaps();
        size_t mref_red_length = mref_length; // reduced length (sites that cannot be quickly imputed)
//        size_t naf0001=0, naf001=0, naf01=0, naf1=0, nafr=0;
        for (size_t mr = mref_start; mr < mref_end; mr++) {
            float af = vcfdata.getAlleleFreqFullRefRegion(mr);
//            if (af <= 0.00001) {
//                naf0001++;
//            } else if (af <= 0.0001) {
//                naf001++;
//            } else if (af <= 0.001) {
//                naf01++;
//            } else if (af <= 0.01) {
//                naf1++;
//            } else
//                nafr++;
            float maf = af > 0.5 ? 1-af : af;
            if (maf*nref < 0.1*m_height) {
                mref_red_length--;
            }
        }
        cout << "      [ " << mref_start << "\t" << mref_end << "\t(" << mref_red_length << "/" << mref_length << " = " << (mref_red_length*100.0/mref_length) << "%) ]" << endl;
                //"\t" << naf0001/(double)mref_length << "\t" << naf001/(double)mref_length << "\t" << naf01/(double)mref_length << "\t" << naf1/(double)mref_length << "\t" << nafr/(double)mref_length << "\t]" << endl;
    }
}
void TargetImp::printSMMatch(const Match &m) {
    cout << "\t" << m.start << "\t" << m.end << "\t(" << (m.end - m.start) << ")\t" <<
                    m.refstart << "\t" << m.refend << "\t(" << (m.refend - m.refstart + 1) << ")" << endl;
}
#endif

// rtree is expected to be NULL on entry!
size_t TargetImp::processHistorySubTree(HistoryTreeNode* &ltree, HistoryTreeNode* &rtree, size_t m, bool curr_hap, bool missing, size_t treedepth, bool perase, vector<Match> &matches) {

    // recursive top-down depth-first processing through tree

    PBWTInterval to0, to1;
    if (!perase) { // need mapping only if the parent will not be erased, otherwise this node will not be extendable anyway
        pbwt.mapInterval(m, ltree->interval, to0, to1);
    }

    bool empty0 = to0.isEmpty();
    bool empty1 = to1.isEmpty();

    // check, how we can extend
    bool extend0 = !empty0 && (missing || !curr_hap || ltree->errors < setmaxerrors);
    bool extend1 = !empty1 && (missing ||  curr_hap || ltree->errors < setmaxerrors);

    // new node with 1-mapping (before processing 0-mapping because we copy the current contents which could be modified below)
    if (extend1) {
        rtree = new HistoryTreeNode(to1, ltree->depth+1, curr_hap||missing ? ltree->errors : ltree->errors+1); // errors plus one if not 1-hap
    }

    if (extend0) { // continue current interval with 0-mapping
        ltree->interval = move(to0);
        ltree->depth++;
        if (curr_hap && !missing)
            ltree->errors++;
    }

    if (!extend0 && !extend1) { // not extendable: empty mapping and no more errors allowed
        // Append potential match if this is a leaf at maximum depth (when the processing of the tree began).
        // There could be other potential matches from other leafs with the same depth as well,
        // or other leafs at the same depth may be extendable.
        // So if this match will be reported as set-maximal must be decided by the caller.
        if (ltree->depth == treedepth && treedepth > 0) // special case at the beginning could report an empty match if there were no matches at all at the first site, which requires the test if treedepth > 0
            matches.emplace_back(m-ltree->depth, m, ltree->interval.getLeft(), ltree->interval.getRight());
    }

    // descend:
    // NOTE: all subtrees have at least the same number of errors, so, for the right extension, if we have not created a new node
    // for rtree, the descending nodes will neither do.

    size_t ldepth = 0, rdepth = 0;

    // descend to left tree
    if (ltree->left)
        ldepth = processHistorySubTree(ltree->left, rtree ? rtree->left : rtree, m, curr_hap, missing, treedepth, !extend0 && !extend1, matches);

    // descend to right tree
    if (ltree->right)
        rdepth = processHistorySubTree(ltree->right, rtree ? rtree->right : rtree, m, curr_hap, missing, treedepth, !extend0 && !extend1, matches);

    size_t depth = ldepth > rdepth ? ldepth : rdepth;

    if (extend0) { // ltree stays alive
        // update depth from ltree
        if (ltree->depth > depth)
            depth = ltree->depth;
        // merge nodes, if mapping is equal
        if (ltree->left && !ltree->right && ltree->left->interval == ltree->interval && ltree->left->errors == ltree->errors) {
            HistoryTreeNode* tmp = ltree;
            ltree=ltree->left;
            delete tmp;
        } else if (!ltree->left && ltree->right && ltree->right->interval == ltree->interval && ltree->right->errors == ltree->errors) {
            HistoryTreeNode* tmp = ltree;
            ltree=ltree->right;
            delete tmp;
        }
    } else { // ltree has to be deleted because we cannnot apply a 0-extension
        delete ltree;
        ltree=NULL; // ltree is a reference to the pointer, so this will set the parent link to NULL
    }

    if (extend1) { // rtree was created
        // update depth from rtree
        if (rtree->depth > depth)
            depth = rtree->depth;
        // merge nodes, if mapping is equal
        if (rtree->left && !rtree->right && rtree->left->interval == rtree->interval && rtree->left->errors == rtree->errors) {
            HistoryTreeNode* tmp = rtree;
            rtree=rtree->left;
            delete tmp;
        } else if (!rtree->left && rtree->right && rtree->right->interval == rtree->interval && rtree->right->errors == rtree->errors) {
            HistoryTreeNode* tmp = rtree;
            rtree=rtree->right;
            delete tmp;
        }
    }

    return depth;
}

void TargetImp::forwardBlock(unsigned block, size_t nsites) {
    // as the imputation function but without actually imputing, only sweeping through the sites and the corresponding matches
    // TODO somehow double code... can this be optimized?

    // all sm matches that cover the current site (i.e. curr >= start && curr <= end-1) are used for imputation
    for (size_t nbunch = 0; nbunch < nsites && mrefs[block] < Mref; nbunch++, mrefs[block]++) {
        if (mrefs[block] == mmaps[block]) { // found common site
            if (beforefirstcommons[block]) {
                // special case: we haven't had a common site yet (first segment), so m needs to stay zero for one more segment
                beforefirstcommons[block] = false;
            } else {
                // set m to this (common) site
                currms[block]++;
            }
            // find first sm match that covers m
            while (idx0s[block] < sm_matches.size() && sm_matches[idx0s[block]].end <= currms[block])
                idx0s[block]++;
        }

        bool missing = currms[block] == nextmisss[block];
        if (missing && mrefs[block] == mmaps[block]) {
            // find next missing
            if (miss_idxs[block] >= misspos.size())
                nextmisss[block] = ~0ull;
            else {
                nextmisss[block] = misspos[miss_idxs[block]];
                miss_idxs[block]++;
            }
        }

        if (mrefs[block] == mmaps[block]) // need to update the mapping for the next common site if we are at a common site
            mmaps[block] = (currms[block]+1) < M ? vcfdata.getSNPIndexToFullRef(currms[block]+1) : Mref;

    }
}

// imputes the next nsites variants according to pre-calculated set-maximal matches and stores the results in imputedTargetBunch and imputedDosageBunch
void TargetImp::imputeBunch(unsigned block, size_t nsites, BooleanVector &imputedTargetBunch, vector<float> &imputedDosageBunch) {
    // iterate over all sites in full reference, keep track of the next common site in both target and reference

    // all sm matches that cover the current site (i.e. curr >= start && curr <= end-1) are used for imputation
    for (size_t nbunch = 0; nbunch < nsites && mrefs[block] < Mref; nbunch++, mrefs[block]++) {

        if (mrefs[block] == mmaps[block]) { // found common site
            if (beforefirstcommons[block]) {
                // special case: we haven't had a common site yet (first segment), so m needs to stay zero for one more segment
                beforefirstcommons[block] = false;
            } else {
                // set m to this (common) site
                currms[block]++;
            }
            // find first sm match that covers m
            while (idx0s[block] < sm_matches.size() && sm_matches[idx0s[block]].end <= currms[block])
                idx0s[block]++;
        }

        bool missing = currms[block] == nextmisss[block];
        if (missing && mrefs[block] == mmaps[block]) {
            // find next missing
            if (miss_idxs[block] >= misspos.size())
                nextmisss[block] = ~0ull;
            else {
                nextmisss[block] = misspos[miss_idxs[block]];
                miss_idxs[block]++;
            }
        }

        // impute only if required
        if (mrefs[block] != mmaps[block] || missing || imputeCalls) {

            if (idx0s[block] >= sm_matches.size() || sm_matches[idx0s[block]].start > currms[block]) { // m points to LAST common site!
                // this site is not covered: As in PBWT we take the reference allele frequency later and count this as a "conflict"
                imputedDosageBunch[nbunch] = DOSAGE_CONFLICT;
            } else {

                // iterate over matches that include the first site
                uint64_t score = 0, sum = 0;
                for (size_t idx = idx0s[block]; idx < sm_matches.size() && sm_matches[idx].start <= currms[block]; idx++) { // m points to LAST common site

                    // single score: lower at the edges of a match, higher in the middle AND the longer the match the higher the score
                    // since we excluded matches of size 1 beforehand, this will always generate a score
                    // (we adapted the score such that m at the beginning of the match will also generate a value)
                    uint64_t single = (currms[block] - sm_matches[idx].start + 1) * (sm_matches[idx].end - currms[block]);

                    // check if quick imputation for this site is possible (if AF is very low)
                    float af = vcfdata.getAlleleFreqFullRefRegion(mrefs[block]);
                    float maf = af > 0.5 ? 1-af : af;
                    size_t sm_height = sm_matches[idx].refend - sm_matches[idx].refstart + 1;
                    // check if due to AF the minor allele in this match cannot contribute enough (less than 0.1 dosage) to the score
                    if (maf * vcfdata.getNReferenceHaps() < 0.01 * sm_height) {
                        sum += sm_height * single;
                        if (af > 0.5) { // if AF is not MAF -> increase score as if every ref would contribute a '1' haplotype
                            score += sm_height * single;
                        } // else, leave score as is (such as every ref would contribute a '0' haplotype)
                    } else {
                        // iterate over all references included in the interval of the match
                        for (size_t ref = sm_matches[idx].refstart; ref <= sm_matches[idx].refend; ref++) {
                            sum += single;
                            size_t absref = pbwt.getSortOrder(sm_matches[idx].end)[ref]; // find the corresponding absolute index in the reference
                            if (vcfdata.getReferenceFullT()[mrefs[block]][vcfdata.getHaploidsRefMap()[absref]]) { // is the corresponding haplotype == 1?
                                score += single;
                            }
                        }
                    }

                }

                // decide on the imputed haplotype according to the score.
                // we can ensure sum > 0 here since matches of size 1 were excluded beforehand and uncovered sites at this point as well
                if (score > sum/2) // only set 1-alleles. if score is even and we have a 50/50, the allele would be 0. That's ok.
                    imputedTargetBunch.setWithPreInit(nbunch, true);
                imputedDosageBunch[nbunch] = (float)(((double) score) / ((double) sum)); // Why double precision first? -> sum and score can get very large!

            }
        }

        if (mrefs[block] == mmaps[block]) { // need to update the mapping for the next common site if we are at a common site
            mmaps[block] = (currms[block]+1) < M ? vcfdata.getSNPIndexToFullRef(currms[block]+1) : Mref;
        }

    }
}
