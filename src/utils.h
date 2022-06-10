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

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <map>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <pwd.h>
#include <sys/inotify.h>
#include <unistd.h>
}

#include "popcntlut.h"
#include "StatusFile.h"

using namespace std;

template<typename T>
string si_binary(T num, unsigned max_int_bits = 1) {
    static const int orders = 9;
    static const char* const symbols[orders] = {
        "B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB"
    };

    T n = num;
    int order;
    for(order = 0; order < orders; order++) {
        if(n <= max_int_bits * 1 * 1024) {
            break;
        }
        n >>= 10;
    }

    if(order < 0)
        order = 0;
    if(order >= orders)
        order = 8;

    num >>= (10*order);
    stringstream ss;
    ss << num << " " << symbols[order];
    return ss.str();
}

inline int popcount64lut(uint64_t p) {
    int res = 0;
//    for (int i = 0; i < 8; i++)
//        res += POPCNT8LUT[(p >> (uint64_t)(8*i)) & 0xff];
    for (int i = 0; i < 4; i++)
        res += POPCNT16LUT[(p >> (uint64_t)(16*i)) & 0xffff];
    return res;
}

template<typename T>
inline T roundToMultiple(T n, T mult) {
    if (n % mult)
        return n + mult - (n % mult);
    else
        return n;
}

template<typename T>
inline T reduceToMultiple(T n, T mult) {
    return n - (n % mult);
}

template<typename T>
inline T divideRounded(T a, T b) {
    return (a + b - 1) / b;
}

// convert NULL terminated C string to a reverse complemented string
// only charactes ATCG are supported, others are kept as they are
inline string reverseComplement(const char *s) {
    string rc(s);
    auto rcit = rc.begin();
    for (int i = rc.size()-1; i >= 0; i--) {
        char c = s[i];
        switch (c) {
        case 'A':
            c = 'T';
            break;
        case 'T':
            c = 'A';
            break;
        case 'C':
            c = 'G';
            break;
        case 'G':
            c = 'C';
            break;
        default:
            break;
        }
        *rcit = c;
        rcit++;
    }
    return rc;
}

inline void addRunlength(int runlength, vector<char>& dest) {
    const int MAXRUNLENGTH = 32767; // 2^31-1
    const int MIN2BYTERUNLENGTH = 128;
    while (runlength > MAXRUNLENGTH) {
        // add maximum run length and following zero byte to encode run lengths exceeding two bytes
        // NOTE: the MSB of the first byte is only an indicator, that the second byte must be considered as well for the runlength
        dest.push_back((char)0xff); // including indicator
        dest.push_back((char)0xff);
        dest.push_back((char)0); // zero byte
        runlength -= MAXRUNLENGTH;
    }
    if (runlength >= MIN2BYTERUNLENGTH) { // runlength needs to be encoded in two bytes
        dest.push_back((char)((runlength & 0x7f) | 0x80)); // lower 7 bits including 2-byte indicator
        dest.push_back((char)(runlength >> 7)); // higher 8 bits
    } else // runlength is encoded in one byte
        dest.push_back((char)runlength);
}

// uses runlength encoding to pack bitssize elements from bits
// dest will be cleared and contains the packed sequence afterwards
// returns true after successful encoding
// NOTE: the encoding stops and returns false, if the encoded sequence will be greater than the source!!
inline bool runlengthEncode(const uint64_t *src, size_t srcsize, vector<char>& dest) {
    dest.clear();
    dest.reserve(srcsize*8+3); // works with any value, but this is probably most efficient for the current implementation; NOTE: srcsize is the number of uint64_t's used for chars
    bool currbit = false;
    int currlength = 0;
    for (size_t i = 0; i < srcsize; i++, src++) {
        uint64_t word = *src;
        // from LSB to MSB
        // NOTE: in order to be able to directly use the builtin_ctzll, word is inverted according to the current bit
        if (currbit)
            word = ~word;
        int rem = 64;
        while (rem > 0) {
            int tr = word == 0 ? rem : __builtin_ctzll(word); // count trailings of current bit
            if (tr < rem) { // add to current runlength and write out
                currlength += tr;
                rem -= tr;
                word >>= tr;
                addRunlength(currlength, dest);
                if (dest.size() >= srcsize*8)
                    return false; // the encoded sequence will be longer than the original, so stop
                currlength = 0;
                currbit = !currbit;
                word = ~word;
            } else { // add remainder to current runlength and continue with next word
                currlength += rem;
                rem = 0;
            }
        }
    }
    // don't forget the last runlength!
    if (currlength)
        addRunlength(currlength, dest);
    if (dest.size() > srcsize*8)
        return false; // the encoded sequence is longer than the original, so return false
    return true;
}

// space at dest must be already reserved for the expected size of the decoded sequence (destsize in number of uint64_t) and pre-initialized with zeros!!!
inline void runlengthDecode(const vector<char> &enc, uint64_t *dest, size_t destsize) {
    bool currbit = true; // start condition: since remlength is 0 it will immediately be switched to false for the first word
    int remlength = 0;
    auto enc_it = enc.cbegin();
//    cout << "start" << endl;
    for (size_t i = 0; i < destsize; i++, dest++) {
        uint64_t *curr = dest;
        int rem = 64;
        while (rem > 0) {
//            cout << "i=" << i << " curr: " << (currbit ? "1" : "0") << " remlength: " << remlength << " rem: " << rem;
            while (remlength == 0) { // get next runlength ("while" is used to handle run lengths of zero correctly)
                if (enc_it == enc.cend()) {
                    StatusFile::addError("Runlength sequence is not long enough for decoding!");
                    exit(EXIT_FAILURE);
                }
                remlength = (*enc_it) & 0x7f; // only the first 7 bits encode the runlength in the first byte
                if ((*enc_it) & 0x80) { // indicator is set, so next byte has to be considered as well
                    enc_it++;
                    if (enc_it == enc.cend()) {
                        StatusFile::addError("Runlength sequence is not long enough for decoding!");
                        exit(EXIT_FAILURE);
                    }
                    remlength |= ((*enc_it) & 0xff) << 7;
                }
                enc_it++;
                currbit = !currbit;
//                cout << "\n load! curr: " << (currbit ? "1" : "0") << " remlength: " << remlength << " rem: " << rem;
            }
            // remlength is > 0 now
            int fill;
            if (remlength >= rem) { // will finish the current word with the current run length
                fill = rem;
                remlength -= rem;
                rem = 0;
            } else { // need to decode more runlength words
                fill = remlength;
                rem -= remlength;
                remlength = 0;
            }
            // fill the next "fill" bits: (only 1's, pre-initialized with 0's)
            if (currbit) {
                // build a vector of "fill" ones
                uint64_t fillvec = fill == 64 ? 0xffffffffffffffffull : ((1ull << fill) - 1ull);
                // shift according to current start position for the fill (64-rem points at the pos after the fill)
                fillvec <<= 64-(rem+fill);
                // apply to word
                *curr |= fillvec;
            }
//            cout << " word: " << hex << setw(16) << setfill('0') << *curr << dec << endl;
        }
    }
    // final check if decoding was correct
    if (remlength || enc_it != enc.cend()) {
        StatusFile::addError("Encoded runlength sequence is longer than excpected!");
        exit(EXIT_FAILURE);
    }
}

inline int setLock(const string& lockdir, const string& lockfile_plain) {

    string lockfile = lockdir + "/" + lockfile_plain;
    mode_t m = umask(0);
    int fd = open(lockfile.c_str(), O_RDWR | O_CREAT, 0666);
    umask(m);
    if(fd >= 0) {
        flock(fd, LOCK_EX);
    } else {
        StatusFile::addError("Could not open lockfile.");
        exit(EXIT_FAILURE);
    }

    // Write some informational file content
    stringstream buf;
    time_t epoch = time(nullptr);
    tm local = *localtime(&epoch);
    buf << "Since: " << put_time(&local, "%c") << endl;

    auto uid = geteuid();
    auto pw = getpwuid(uid);
    buf << "By: " << pw->pw_name << " (" << pw->pw_gecos << ")" << endl;

    buf << "Pid: " << getpid() << endl;
    char *cmdline = new char[4096]; // arbitrary
    int cmdline_fd = open("/proc/self/cmdline", O_RDONLY);
    if(cmdline_fd != -1) {
        cmdline[0] = 0;
        cmdline[4095] = 0;
        if(read(cmdline_fd, cmdline, 4096) == -1) {
            StatusFile::addError("Could not read process details from /proc/self");
            exit(EXIT_FAILURE);
        }
        buf << "Cmd: " << cmdline << endl;
    }
    delete[] cmdline;

    if(write(fd, buf.str().c_str(), buf.str().length()) == -1) {
        StatusFile::addError("Could not write to lockfile.");
        exit(EXIT_FAILURE);
    }

    return fd;

}

inline void releaseLock(int fd) {
    int ftr = ftruncate(fd, 0); // clear the contents of the lockfile
    if (ftr != 0)
        StatusFile::addWarning("Release lock returned status " + to_string(ftr) + ".");
    close(fd);
}

#endif // UTILS_H
