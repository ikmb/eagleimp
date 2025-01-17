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

#ifndef MYMALLOC_H_
#define MYMALLOC_H_

#include <cstdlib>
#include <map>
#include <atomic>
#include <mutex>
#include <string>
#include <utility>

using namespace std;

class MyMalloc {
public:

    static void *malloc(size_t size, const string &id) {
        return mm.mymalloc(size, id);
    }

    static void *calloc(size_t size, size_t elt_size, const string &id) {
        return mm.mycalloc(size, elt_size, id);
    }

    static void *realloc(void *p, size_t size, const string &id) {
        return mm.myrealloc(p, size, id);
    }

    static void free(void *p) {
        mm.myfree(p);
    }

    // don't care about mutual exclusion here
    static size_t getTotalAlloced() { return mm.mem_alloced; }
    static size_t getTotalFreed() { return mm.mem_freed; }
    static size_t getCurrentAlloced() { return mm.mem_alloced - mm.mem_freed; }
    static size_t getMaxAlloced() { return mm.mem_maxalloced; }
    static const map<size_t,pair<size_t,string>> &getMap() { return mm.mem_map; };

    static void printSummary(const string &note) {
        cerr << "\nMemory " << note << ":" << endl;
        cerr << "allocated:    " << getTotalAlloced() << " (" << MyMalloc::getTotalAlloced()/(1024*1024) << " MiB)" << endl;
        cerr << "freed:        " << getTotalFreed() << " (" << MyMalloc::getTotalFreed()/(1024*1024) << " MiB)" << endl;
        cerr << "not freed:    " << getCurrentAlloced() << " (" << MyMalloc::getCurrentAlloced()/(1024*1024) << " MiB)" << endl;
        cerr << "max. alloced: " << getMaxAlloced() << " (" << MyMalloc::getMaxAlloced()/(1024*1024) << " MiB)" << endl;
        const auto &memmap = mm.mem_map;
        cerr << "Remaining map entries: " << memmap.size() << endl;
        if (memmap.size()) {
            int count_refdata = 0;
            int count_refdataT = 0;
            int count_refdataOverlap = 0;
            int count_refdataTOverlap = 0;
            int count_hapdata = 0;
            int count_tgt_gt = 0;
            int count_else = 0;
            size_t size_refdata = 0;
            size_t size_refdataT = 0;
            size_t size_refdataOverlap = 0;
            size_t size_refdataTOverlap = 0;
            size_t size_hapdata = 0;
            size_t size_tgt_gt = 0;
            size_t size_else = 0;
            for (const auto &m : memmap) {
                if (m.second.second.find("refdataTOverlap") != string::npos) {
                    count_refdataTOverlap++;
                    size_refdataTOverlap += m.second.first;
                } else if (m.second.second.find("refdataOverlap") != string::npos) {
                    count_refdataOverlap++;
                    size_refdataOverlap += m.second.first;
                } else if (m.second.second.find("refdataT") != string::npos) {
                    count_refdataT++;
                    size_refdataT += m.second.first;
                } else if (m.second.second.find("refdata") != string::npos) {
                    count_refdata++;
                    size_refdata += m.second.first;
                } else if (m.second.second.find("hapdata") != string::npos) {
                    count_hapdata++;
                    size_hapdata += m.second.first;
                } else if (m.second.second.find("tgt_gt") != string::npos) {
                    count_tgt_gt++;
                    size_tgt_gt += m.second.first;
                } else {
                    count_else++;
                    size_else += m.second.first;
                }
            }
            if (count_refdata)
                cerr << "  refdata:         " << count_refdata << "\t" << size_refdata / 1024 << " kiB" << endl;
            if (count_refdataT)
                cerr << "  refdataT:        " << count_refdataT << "\t" << size_refdataT / 1024 << " kiB" << endl;
            if (count_refdataOverlap)
                cerr << "  refdataOverlap:  " << count_refdataOverlap << "\t" << size_refdataOverlap / 1024 << " kiB" << endl;
            if (count_refdataTOverlap)
                cerr << "  refdataTOverlap: " << count_refdataTOverlap << "\t" << size_refdataTOverlap / 1024 << " kiB" << endl;
            if (count_hapdata)
                cerr << "  hapdata:         " << count_hapdata << "\t" << size_hapdata / 1024 << " kiB" << endl;
            if (count_tgt_gt)
                cerr << "  tgt_gt:          " << count_tgt_gt << "\t" << size_tgt_gt / 1024 << " kiB" << endl;
            if (count_else)
                cerr << "  other:           " << count_else << "\t" << size_else / 1024 << " kiB" << endl;
        }
    }

    static void dumpMemMap(ofstream &ofs) {
        const auto &memmap = mm.mem_map;
        if (memmap.size()) {
            ofs << "not freed:    " << MyMalloc::getCurrentAlloced() << " (" << MyMalloc::getCurrentAlloced()/1024 << " kiB)" << endl;
            ofs << "max. alloced: " << MyMalloc::getMaxAlloced() << " (" << MyMalloc::getMaxAlloced()/1024 << " kiB)" << endl;
            ofs << "Remaining map entries (size in kiB): " << memmap.size() << endl;
            for (const auto &m : memmap)
                ofs << m.second.second << "\t" << m.second.first/1024 << endl;
        }
    }


private:
    MyMalloc(){}

    static MyMalloc mm;

    void *mymalloc(size_t size, const string &id __attribute__((unused))) {
        void *p = std::malloc(size);
        add(p, size, id);
        return p;
    }

    void *mycalloc(size_t size, size_t elt_size, const string &id __attribute__((unused))) {
        void *p = std::calloc(size, elt_size);
        add(p, size, id);
        return p;
    }

    void *myrealloc(void *p, size_t size, const string &id __attribute__((unused))) {
        void *oldp __attribute__((unused)) = p;
        void *np = std::realloc(p, size);
        if (np) { // re-allocation successful
            remove(oldp);
            add(np, size, id);
        }
        return np;
    }

    void myfree(void *p) {
        remove(p);
        std::free(p);
    }

    inline void add(void *p, size_t size, const string &id) {
        if (p) {
            mux.lock();
            mem_alloced += size;
            if (mem_alloced-mem_freed > mem_maxalloced) {
//                cerr << "MyMalloc: reached max! added:\t" << size << "\t" << id << "\tin use: " << mem_alloced - mem_freed << endl;
                mem_maxalloced = mem_alloced - mem_freed;
            }
//            if (mem_map.find((size_t)p) != mem_map.end())
//                cerr << "MyMalloc: Address already in use! Double allocation?! ID: " << id << endl;
            mem_map[(size_t)p] = make_pair(size,id);
            mux.unlock();
        }
    }

    inline void remove(void *p) {
        mux.lock();
        size_t sizefreed = mem_map[(size_t)p].first;
        mem_freed += sizefreed;
//        if (sizefreed == 0)
//            cerr << "MyMalloc: freed size is zero. Did you allocate with MyMalloc?" << endl;
        mem_map.erase((size_t)p);
        mux.unlock();
    }

    mutex mux;
    map<size_t,pair<size_t,string>> mem_map;
    size_t mem_alloced = 0;
    size_t mem_freed = 0;
    size_t mem_maxalloced = 0;

};
#endif /* MYMALLOC_H_ */
