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

#ifndef STOPWATCH_H_
#define STOPWATCH_H_

#include <string>
#include <sstream>
#include <vector>
#include <array>

extern "C" {
    #include <sys/resource.h>
}

/**********************/
/* NOT THREADSAFE !!! */
/**********************/

class Stopwatch {
private:
    bool active = true;

    class StopwatchGlobal {
    public:

        // no copies!
        StopwatchGlobal(const Stopwatch&) = delete;
        StopwatchGlobal &operator=(const Stopwatch&) = delete;

        struct res_info {
            rusage r_laststart = {{0,0},{0,0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0}};
            rusage r_elapsed = {{0,0},{0,0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0}};
            timeval t_laststart = {0,0};
            timeval t_elapsed = {0,0};
            std::string id;
            bool active;
        };

        static StopwatchGlobal &getInstance() {
            return instance;
        }

        // start a timer with the given identifier
        void start(const char *ident);
        // stop last timer
        void stop();
        // dump all timers
        void dump(std::ostream &out);

    private:
        template <class T> class Tree {
        public:
            Tree(T content_) : content(content_) {}
            void addChild(T child_content) { children.push_back(Tree(child_content)); }
            std::vector<Tree<T>>& getChildren() { return children; }
            T& getContent() { return content; }
            void setContent(T& content_) { content = content_; }
        private:
            T content;
            std::vector<Tree<T>> children;
        };

        static res_info createNewResInfo(const char* id);

        StopwatchGlobal();

        void dump(std::ostream &out, Tree<res_info>& subtree, int depth, std::vector<bool>& structinfo, std::array<int,10>& colwidths);

        static void print_resinfo(std::ostream &out, const res_info &info, int depth, std::vector<bool>& structinfo, std::array<int,10>& colwidths);

        Tree<res_info> resources;
        std::vector<Tree<res_info>*> lastactivated;
        int maxnamelen;

        static StopwatchGlobal instance;

        static constexpr const char* ROOTNAME = "**TOTAL**";
        static constexpr const size_t MINCOLWIDTH = 8;
        static constexpr const size_t MINCOLSEP = 1;
    }; // END private class StopwatchGlobal

public:
    // no default constructor
    Stopwatch() = delete;
    // no copies!
    Stopwatch(const Stopwatch&) = delete;
    Stopwatch &operator=(const Stopwatch&) = delete;

    // explicit constructor with name, starts immediately
    explicit Stopwatch(const char* ident) {
        StopwatchGlobal::getInstance().start(ident);
    }

    // stop this stopwatch manually
    void stop() {
        if (!active)
            std::cerr << "WARNING! Stopwatch stopped twice!" << std::endl;
        else {
            StopwatchGlobal::getInstance().stop();
            active = false;
        }
    }

    // destructor stops this stopwatch (if not done already manually)
    ~Stopwatch() { if (active) stop(); }

    static void dump(std::ostream &out) {
        StopwatchGlobal::getInstance().dump(out);
    }

};

#endif /* STOPWATCH_H_ */
