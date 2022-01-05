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

#include <iostream>
#include <iomanip>
#include <string>
#include <sys/resource.h>
#include <sys/time.h>
#include <algorithm>

#include "Stopwatch.h"

using namespace std;

Stopwatch::StopwatchGlobal Stopwatch::StopwatchGlobal::instance;

Stopwatch::StopwatchGlobal::StopwatchGlobal() : resources(createNewResInfo(ROOTNAME)), maxnamelen(string(ROOTNAME).size()) {
    // add root
    lastactivated.push_back(&resources);
}

/**
 * Calculates the difference between time values (x -= y).
 * @param x
 * @param y (const)
 */
static inline void tv_diff(timeval &x, const timeval &y) {
    x.tv_sec -= y.tv_sec;
    x.tv_usec -= y.tv_usec;
    if (x.tv_usec < 0) {
        --x.tv_sec;
        x.tv_usec += 1000000;
    }
}

/**
 * Accumulates time values (x += y).
 * @param x
 * @param y (const)
 */
static inline void tv_acc(timeval &x, const timeval &y) {
    x.tv_sec += y.tv_sec;
    x.tv_usec += y.tv_usec;
    if (x.tv_usec >= 1000000) {
        ++x.tv_sec;
        x.tv_usec -= 1000000;
    }
}

/**
 * Calculates the difference between rusage values (x -= y). (Only the values required here.)
 * @param x
 * @param y (const)
 */
static inline void rusage_diff(rusage &x, const rusage &y) {
    tv_diff(x.ru_utime, y.ru_utime);
    tv_diff(x.ru_stime, y.ru_stime);
    x.ru_inblock -= y.ru_inblock;
    x.ru_oublock -= y.ru_oublock;
    x.ru_minflt  -= y.ru_minflt;
    x.ru_majflt  -= y.ru_majflt;
    x.ru_nvcsw   -= y.ru_nvcsw;
    x.ru_nivcsw  -= y.ru_nivcsw;
}

/**
 * Accumulates rusage values (x += y). (Only the values required here.)
 * @param x
 * @param y (const)
 */
static inline void rusage_acc(rusage &x, const rusage &y) {
    tv_acc(x.ru_utime, y.ru_utime);
    tv_acc(x.ru_stime, y.ru_stime);
    x.ru_inblock += y.ru_inblock;
    x.ru_oublock += y.ru_oublock;
    x.ru_minflt  += y.ru_minflt;
    x.ru_majflt  += y.ru_majflt;
    x.ru_nvcsw   += y.ru_nvcsw;
    x.ru_nivcsw  += y.ru_nivcsw;
}

static inline size_t stringsize(unsigned val) {
    stringstream s;
    s << val;
    return s.str().size();
}

/* static */ Stopwatch::StopwatchGlobal::res_info Stopwatch::StopwatchGlobal::createNewResInfo(const char* id) {
    res_info info;
    getrusage(RUSAGE_SELF, &(info.r_laststart));
    gettimeofday(&(info.t_laststart), NULL);
    info.id = string(id);
    info.active = true;
    return info;
}

void Stopwatch::StopwatchGlobal::start(const char *ident) {

    auto res_iter = lastactivated.back()->getChildren().begin();

    for (; res_iter != lastactivated.back()->getChildren().end(); res_iter++) {
        if (res_iter->getContent().id.compare(ident) == 0) { // found resource with same identifier -> restart
            res_info &info = res_iter->getContent();
            if (!info.active) {
                getrusage(RUSAGE_SELF, &(info.r_laststart));
                gettimeofday(&(info.t_laststart), NULL);
                info.active = true;
                lastactivated.push_back(&(*res_iter));
        //        std::cout << "RESTART: ";
        //        print_resinfo(std::cout, ident, info);
            } // else nothing to do, is already active
            break;
        }
    }

    if(res_iter == lastactivated.back()->getChildren().end()) { // new timer required
        maxnamelen = max(maxnamelen, (int)(string(ident).size() + lastactivated.size()));
        lastactivated.back()->addChild(createNewResInfo(ident));
        lastactivated.push_back(&(lastactivated.back()->getChildren().back()));
//        std::cout << "START: ";
//        print_resinfo(std::cout, ident, info);
    }
}

void Stopwatch::StopwatchGlobal::stop() {

    if(lastactivated.size() == 0)
        throw runtime_error("ERROR! Attempt no timers to stop!");

    res_info &info = lastactivated.back()->getContent();

    if (info.active) { // only if the timer is active, but this has to be the case!
        rusage u;
        getrusage(RUSAGE_SELF, &u);

        timeval tv;
        gettimeofday(&tv, NULL);

        tv_diff(tv, info.t_laststart);
        rusage_diff(u, info.r_laststart);

        tv_acc(info.t_elapsed, tv);
        rusage_acc(info.r_elapsed, u);

        info.active = false;

        lastactivated.resize(lastactivated.size()-1); // remove last element
    } else // not active -> error!
        throw runtime_error("ERROR: inactive timer claimed to be active!");
//    std::cout << "STOP: ";
//    print_resinfo(std::cout, ident, info);
}

void Stopwatch::StopwatchGlobal::dump(std::ostream &out) {
    // first stop all running timers (all but the root timer should have already been stopped...)
    while (lastactivated.size() > 0)
        stop();

    // determine column widths
    // name
    array<int,10> colwidths;
    colwidths[0] = maxnamelen + MINCOLSEP;
    // for all other columns, the maximum values (with thus, the largest sizes) are in the root timer
    res_info& root = resources.getContent();
    // for the time columns, the minimum column width is reduced by 3 precision digits and the decimal point
    colwidths[1] = max(MINCOLWIDTH-4, stringsize(root.t_elapsed.tv_sec) + MINCOLSEP);
    colwidths[2] = max(MINCOLWIDTH-4, stringsize(root.r_elapsed.ru_utime.tv_sec) + MINCOLSEP);
    colwidths[3] = max(MINCOLWIDTH-4, stringsize(root.r_elapsed.ru_stime.tv_sec) + MINCOLSEP);
    // the minimum column width for other columns is not reduced
    colwidths[4] = max(MINCOLWIDTH, stringsize(root.r_elapsed.ru_inblock) + MINCOLSEP);
    colwidths[5] = max(MINCOLWIDTH, stringsize(root.r_elapsed.ru_oublock) + MINCOLSEP);
    colwidths[6] = max(MINCOLWIDTH, stringsize(root.r_elapsed.ru_minflt) + MINCOLSEP);
    colwidths[7] = max(MINCOLWIDTH, stringsize(root.r_elapsed.ru_majflt) + MINCOLSEP);
    colwidths[8] = max(MINCOLWIDTH, stringsize(root.r_elapsed.ru_nvcsw) + MINCOLSEP);
    colwidths[9] = max(MINCOLWIDTH, stringsize(root.r_elapsed.ru_nivcsw) + MINCOLSEP);

    // print header
    // left bound
    out << "TIMER";
    for (int i = 0; i < colwidths[0]-5; i++)
        out << " ";
    // right bound
    for (int i = 0; i < colwidths[1]-3; i++)
        out << " ";
    out << "ELAPSED";
    for (int i = 0; i < colwidths[2]-3; i++)
        out << " ";
    out << "CPUTIME";
    for (int i = 0; i < colwidths[3]-3; i++)
        out << " ";
    out << "SYSTIME";
    for (int i = 0; i < colwidths[4]-6; i++)
        out << " ";
    out << "IOREAD";
    for (int i = 0; i < colwidths[5]-7; i++)
        out << " ";
    out << "IOWRITE";
    for (int i = 0; i < colwidths[6]-6; i++)
        out << " ";
    out << "MINFLT";
    for (int i = 0; i < colwidths[7]-6; i++)
        out << " ";
    out << "MAJFLT";
    for (int i = 0; i < colwidths[8]-5; i++)
        out << " ";
    out << "VCONT";
    for (int i = 0; i < colwidths[9]-7; i++)
        out << " ";
    out << "INVCONT" << endl;

    vector<bool> structinfo; // empty structural information for root

    // recursively dump all trees, start with the root
    dump(out, resources, 0, structinfo, colwidths);

//    for (int c : colwidths)
//        cout << " " << c;
//    cout << endl;
}

void Stopwatch::StopwatchGlobal::dump(std::ostream &out, Tree<res_info>& subtree, int depth, vector<bool>& structinfo, array<int,10>& colwidths) {

    print_resinfo(out, subtree.getContent(), depth, structinfo, colwidths);

    // step down in the tree
    auto& children = subtree.getChildren();

    if (!children.empty())
        structinfo.resize(depth+1); // size according to current depth
    for (size_t i = 0; i < children.size(); i++) {
        if (i == children.size()-1) // last child
            structinfo[depth] = false;
        else
            structinfo[depth] = true;
        dump(out, children[i], depth+1, structinfo, colwidths);
    }
}

/*static*/ void Stopwatch::StopwatchGlobal::print_resinfo(ostream &out, const res_info &info, int depth, vector<bool>& structinfo, array<int,10>& colwidths) {

    // print leading spaces/lines according to depth and structural information
    for (int i = 0; i < depth-1; i++) {
        if (structinfo[i])
            out << "│";
        else
            out << " ";
    }
    if (depth > 0) {
        if (structinfo[depth-1])
            out << "├";
        else
            out << "└";
    }

    // timer name
    out << info.id;
    for (int i = 0; i < colwidths[0] - ((int)info.id.size() + depth); i++)
        out << " ";

    // elapsed time
    out << setw(colwidths[1]) << info.t_elapsed.tv_sec << "." << setfill('0') << setw(3) << (info.t_elapsed.tv_usec / 1000) << setfill(' ');

    // CPU (user) time
    out << setw(colwidths[2]) << info.r_elapsed.ru_utime.tv_sec << "." << setfill('0') << setw(3) << (info.r_elapsed.ru_utime.tv_usec / 1000) << setfill(' ');

    // System (IO) time
    out << setw(colwidths[3]) << info.r_elapsed.ru_stime.tv_sec << "." << setfill('0') << setw(3) << (info.r_elapsed.ru_stime.tv_usec / 1000) << setfill(' ');

    // I/O reads
    out << setw(colwidths[4]) << (unsigned) (info.r_elapsed.ru_inblock);

    // I/O writes
    out << setw(colwidths[5]) << (unsigned) (info.r_elapsed.ru_oublock);

    // minor page faults
    out << setw(colwidths[6]) << (unsigned) (info.r_elapsed.ru_minflt);

    // major page faults
    out << setw(colwidths[7]) << (unsigned) (info.r_elapsed.ru_majflt);

    // voluntary context switches
    out << setw(colwidths[8]) << (unsigned) (info.r_elapsed.ru_nvcsw);

    // involuntary context switches
    out << setw(colwidths[9]) << (unsigned) (info.r_elapsed.ru_nivcsw);

    out << endl;

}
