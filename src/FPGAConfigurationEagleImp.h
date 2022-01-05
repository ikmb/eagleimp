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

#ifndef FPGACONFIGURATIONEAGLEIMP_H_
#define FPGACONFIGURATIONEAGLEIMP_H_

#include <tuple>
#include <stdexcept>
#include <sstream>

#include "hybridsys/Hybridsys.h"
#include "hybridsys/FPGA.h"

using namespace std;

class FPGAConfigurationEagleImp {
public:
    enum ApplicationType {
        APP_INVALID = 0,
        APP_EAGLEIMP  = 4
    };

    FPGAConfigurationEagleImp(){}

    FPGAConfigurationEagleImp(const void *buffer) {
        if (buffer != NULL)
            parse(buffer);
    }

    void parse(const void *buffer);

    unsigned getAppID() const {
        return appID;
    }
    tuple<int,int,int> getVersion() const {
        return make_tuple(versionMajor,versionMinor,versionRevision);
    }
    uint64_t getCrefFrequency() const {
        return crefFrequency;
    }
    uint64_t getPBWTFrequency() const {
        return pbwtFrequency;
    }
    uint64_t getAvailableRAM() const {
        return availableRAM;
    }
    unsigned getNumPipelines() const {
        return numPipelines;
    }
    size_t getMaxSites() const {
        return maxSites;
    }
    size_t getMaxHaps() const {
        return maxHaps;
    }
    size_t getMaxK() const {
        return maxK;
    }

    bool operator==(const FPGAConfigurationEagleImp& other) const;

    static constexpr int major_required = 0x1;
    static constexpr int minor_required = 0x5;

private:
    unsigned appID = APP_INVALID;
    int versionMajor = 0;
    int versionMinor = 0;
    int versionRevision = 0;
    uint64_t crefFrequency = 1; // don't use 0 to prevent division by zero
    uint64_t pbwtFrequency = 1; // don't use 0 to prevent division by zero
    uint64_t availableRAM = 0;
    unsigned numPipelines = 1;
    size_t maxSites = 0;
    size_t maxHaps = 0;
    size_t maxK = 0;

};

ostream& operator<<(ostream& out, const FPGAConfigurationEagleImp &c);

vector<FPGAConfigurationEagleImp> verifyConfiguration(hybridsys::Hybridsys& hs, bool &usefpga, bool &usegpu, bool debug);

#endif /* FPGACONFIGURATIONEAGLEIMP_H_ */
