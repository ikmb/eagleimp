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

#define DISABLE_MYMALLOC

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

#ifndef DISABLE_MYMALLOC
    // don't care about mutual exclusion here
    static size_t getAlloced() { return mm.mem_alloced; }
    static size_t getFreed() { return mm.mem_freed; }
    static size_t getNotFreed() { return mm.mem_alloced - mm.mem_freed; }
    static size_t getMaxAlloced() { return mm.mem_maxalloced; }
    static const map<size_t,pair<size_t,string>> &getMap() { return mm.mem_map; };

    static void print(ofstream &ofs, const string &note) {
        cerr << "\nMemory " << note << ":" << endl;
        cerr << "allocated:    " << getAlloced() << " (" << MyMalloc::getAlloced()/(1024*1024) << " MiB)" << endl;
        cerr << "freed:        " << getFreed() << " (" << MyMalloc::getFreed()/(1024*1024) << " MiB)" << endl;
        cerr << "not freed:    " << getNotFreed() << " (" << MyMalloc::getNotFreed()/(1024*1024) << " MiB)" << endl;
        cerr << "max. alloced: " << getMaxAlloced() << " (" << MyMalloc::getMaxAlloced()/(1024*1024) << " MiB)" << endl;
        const auto &memmap = mm.mem_map;
        cerr << "Remaining map entries (size in MiB): " << memmap.size() << endl;
        if (memmap.size()) {
            auto mapit = memmap.cbegin();
            for (int i=0; i < 20 && mapit != memmap.cend(); i++, mapit++) {
                cerr << " " << i << ":\t" << mapit->second.second << "\t" << mapit->second.first/(1024*1024) << endl;
            }
            if (mapit != memmap.cend())
                cerr << " ..." << endl;
            ofs << "not freed:    " << MyMalloc::getNotFreed() << " (" << MyMalloc::getNotFreed()/1024 << " kiB)" << endl;
            ofs << "max. alloced: " << MyMalloc::getMaxAlloced() << " (" << MyMalloc::getMaxAlloced()/1024 << " kiB)" << endl;
            ofs << "Remaining map entries (size in kiB): " << memmap.size() << endl;
            for (const auto &m : memmap)
                ofs << m.second.second << "\t" << m.second.first/1024 << endl;
        }
    }

#endif

private:
    MyMalloc(){}

    static MyMalloc mm;

    void *mymalloc(size_t size, const string &id __attribute__((unused))) {
        void *p = std::malloc(size);
#ifndef DISABLE_MYMALLOC
        add(p, size, id);
#endif
        return p;
    }

    void *mycalloc(size_t size, size_t elt_size, const string &id __attribute__((unused))) {
        void *p = std::calloc(size, elt_size);
#ifndef DISABLE_MYMALLOC
        add(p, size, id);
#endif
        return p;
    }

    void *myrealloc(void *p, size_t size, const string &id __attribute__((unused))) {
        void *oldp __attribute__((unused)) = p;
        void *np = std::realloc(p, size);
#ifndef DISABLE_MYMALLOC
        if (np) { // re-allocation successful
            remove(oldp);
            add(np, size, id);
        }
#endif
        return np;
    }

    void myfree(void *p) {
#ifndef DISABLE_MYMALLOC
        remove(p);
#endif
        std::free(p);
    }

#ifndef DISABLE_MYMALLOC
    inline void add(void *p, size_t size, const string &id) {
        if (p) {
            mux.lock();
            if (mem_map.find((size_t)p) != mem_map.end())
                cerr << "MyMalloc: Address already in use! Double allocation?! ID: " << id << endl;
            mem_alloced += size;
            mem_map[(size_t)p] = make_pair(size,id);
            if (mem_alloced-mem_freed > mem_maxalloced) {
//                cerr << "MyMalloc: reached max! added:\t" << size << "\t" << id << "\tin use: " << mem_alloced - mem_freed << endl;
                mem_maxalloced = mem_alloced - mem_freed;
            }
            mux.unlock();
        }
    }

    inline void remove(void *p) {
        mux.lock();
        size_t sizefreed = mem_map[(size_t)p].first;
        mem_freed += sizefreed;
        if (sizefreed == 0)
            cerr << "MyMalloc: freed size is zero. Did you allocate with MyMalloc?" << endl;
        mem_map.erase((size_t)p);
        mux.unlock();
    }

    mutex mux;
    map<size_t,pair<size_t,string>> mem_map;
    size_t mem_alloced = 0;
    size_t mem_freed = 0;
    size_t mem_maxalloced = 0;
#endif

};
#endif /* MYMALLOC_H_ */
