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

#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <tbb/concurrent_queue.h>
#include <thread>
#include <vector>
#include <functional>
#include <atomic>

template<class Tin, class Tout>
class ThreadPool
{
public:
    using inqueue_type = tbb::concurrent_bounded_queue<Tin>;
    using outqueue_type = tbb::concurrent_bounded_queue<Tout>;

    using threadfunc_type = std::function<void(inqueue_type& inqueue, outqueue_type& outqueue, std::atomic_bool& terminationRequest, int threadIndex)>;

    ThreadPool(unsigned numThreads_,
               inqueue_type& inqueue_,
               outqueue_type& outqueue_,
               threadfunc_type threadFunc_
               ) :
        numThreads(numThreads_),
        inqueue(inqueue_),
        outqueue(outqueue_),
        threadFunc(threadFunc_),
        terminationRequest(false)
    {}

    ~ThreadPool() {
        setTerminationRequest(true);
        waitForTermination();
    }

    void run() {
        for(unsigned i = 0; i < numThreads; i++) {
            threads.emplace_back(threadFunc, std::ref(this->inqueue), std::ref(this->outqueue), std::ref(this->terminationRequest), i);
        }
    }

    // set the indicator for a termination request.
    // if inqueue_abort is set, abort the inqueue as well.
    void setTerminationRequest(bool inqueue_abort) {
        // actually cancel threads
        terminationRequest = true;
        if (inqueue_abort)
            inqueue.abort();
    }

    void waitForTermination() {
        for(std::thread &t : threads)
            t.join();
        threads.clear();
    }

private:
    std::vector<std::thread> threads;
    unsigned numThreads;
    inqueue_type& inqueue;
    outqueue_type& outqueue;
    threadfunc_type threadFunc;
    std::atomic_bool terminationRequest;

};

#endif /* THREADPOOL_H_ */
