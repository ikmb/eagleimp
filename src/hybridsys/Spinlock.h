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

#ifndef SPINLOCK_H
#define SPINLOCK_H

#include <atomic>

/**
 * @brief A lightweight Spinlock
 *
 * This Spinlock supports the concepts BasicLockable and Lockable and is
 * therefore usable with std::lock_guard. Use this if you are only locking
 * for a very short time. Locking is very efficient (contrary to std::mutex),
 * and waiting is implemented with busy waiting (also contrary to std::mutex).
 *
 * @todo Consider implement this using test-and-test-and-set protocol to reduce
 * performance cost in high contention situations
 */
class Spinlock {
private:
    std::atomic_flag flag = ATOMIC_FLAG_INIT;
public:
    void lock() noexcept {
        while(flag.test_and_set(std::memory_order_acquire)) ;
    }

    void unlock() noexcept {
        flag.clear(std::memory_order_release);
    }

    bool try_lock() noexcept {
        return !(flag.test_and_set(std::memory_order_acquire));
    }
};

#endif // SPINLOCK_H
