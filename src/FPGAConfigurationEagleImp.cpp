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

#include <stdint.h>

#include "hybridsys/FPGA.h"

#include "StatusFile.h"

#include "FPGAConfigurationEagleImp.h"

void FPGAConfigurationEagleImp::parse(const void *buffer) {
    const uint16_t *info = reinterpret_cast<const uint16_t*>(buffer);

    appID            = info[0] >> 8;
    versionMajor     = info[0] & 0xff;
    versionMinor     = info[1] >> 8;
    versionRevision  = info[1] & 0xff;
    crefFrequency    = (((uint32_t)info[3]) << 16) | ((uint32_t) info[2]);
    pbwtFrequency    = (((uint32_t)info[5]) << 16) | ((uint32_t) info[4]);
    availableRAM     = ((((uint64_t)info[7]) << 16) | ((uint64_t)info[6])) << 20; // convert MiB to B
    numPipelines     = info[8] & 0xff;
    maxSites         = info[9];
    maxHaps          = info[10];
    maxK             = info[11];
    // other fields are reserved and
    // info[15] includes the status byte
    // need to multiply the maxX values by 1024:
    maxSites <<= 10;
    maxHaps <<= 10;
    maxK <<= 10;
}

bool FPGAConfigurationEagleImp::operator==(const FPGAConfigurationEagleImp &other) const {
    return
        (appID            == other.appID           ) &&
        (versionMajor     == other.versionMajor    ) &&
        (versionMinor     == other.versionMinor    ) &&
        (versionRevision  == other.versionRevision ) &&
        (crefFrequency    == other.crefFrequency   ) &&
        (pbwtFrequency    == other.pbwtFrequency   ) &&
        (availableRAM     == other.availableRAM    ) &&
        (numPipelines     == other.numPipelines    ) &&
        (maxSites         == other.maxSites        ) &&
        (maxHaps          == other.maxHaps         ) &&
        (maxK             == other.maxK            );
}

std::ostream& operator<<(std::ostream& out, const FPGAConfigurationEagleImp &c) {
    out << "\tApplication type:             ";
    switch(c.getAppID()) {
    case FPGAConfigurationEagleImp::APP_EAGLEIMP: out << "EagleImp"; break;
    default: out << "unknown"; break;
    }
    out << std::endl;
    out << "\tVersion:                      " << std::get<0>(c.getVersion()) << "." << std::get<1>(c.getVersion()) << "." << std::get<2>(c.getVersion()) << std::endl;
    out << "\tFrequency of cref pipeline:   " << c.getCrefFrequency()/1000000 << " MHz" << std::endl;
    out << "\tFrequency of PBWT pipeline:   " << c.getPBWTFrequency()/1000000 << " MHz" << std::endl;
    out << "\tAvailable RAM:                " << c.getAvailableRAM()/1048576 << " MiB" << std::endl;
    out << "\tNumber of pipelines:          " << c.getNumPipelines() << std::endl;
    out << "\tmaximum no. of sites:         " << c.getMaxSites() << std::endl;
    out << "\tmaximum no. of haps:          " << c.getMaxHaps() << std::endl;
    out << "\tmaximum K:                    " << c.getMaxK() << std::endl;
    return out;
}

vector<FPGAConfigurationEagleImp> verifyConfiguration(hybridsys::Hybridsys& hs, bool &usefpga, bool &usegpu, bool debug) {

    if(hs.getFPGAs().size() > 1) {
        StatusFile::addError("Using more than one FPGA for one target file has not yet been implemented. Sorry...");
        exit(EXIT_FAILURE);
    }

    if (hs.getGPUs().size() > 0 && hs.getFPGAs().size() > 0) {
        StatusFile::addError("Using GPUs and FPGAs together has not yet been implemented. Sorry...");
        exit(EXIT_FAILURE);
    }

    usefpga = hs.getFPGAs().size() > 0;
    usegpu  = hs.getGPUs().size() > 0;

    // get all FPGA configuration blocks
    vector<FPGAConfigurationEagleImp> configs;
    configs.resize(hs.getFPGAs().size());

    if (usefpga) {

        for(unsigned i = 0; i < configs.size(); i++) {
            hybridsys::FPGA &f = hs.getFPGA(i);
            hs.getFPGAs()[i].createThreadHandle();
            configs[i].parse(f.getConfigurationData());
            hs.getFPGAs()[i].destroyThreadHandle();
        }

        // verify the FPGA is configured with the correct application
        if (configs[0].getAppID() != FPGAConfigurationEagleImp::APP_EAGLEIMP) {
            StatusFile::addError("The selected FPGA is not configured with the EagleImp application.");
            exit(EXIT_FAILURE);
        }

        // check version
        for (auto &cfg: configs) {
            if (get<0>(cfg.getVersion()) != FPGAConfigurationEagleImp::major_required) {
                stringstream ss;
                ss << "Found FPGA version " << get<0>(cfg.getVersion()) << "." << get<1>(cfg.getVersion()) << "rev" << get<2>(cfg.getVersion())
                   << " which is incompatible with this application. Requiring major version " << FPGAConfigurationEagleImp::major_required << " (at least v"
                   << FPGAConfigurationEagleImp::major_required << "." << FPGAConfigurationEagleImp::minor_required << ").";
                StatusFile::addError(ss.str());
                exit(EXIT_FAILURE);
            }
            if (get<1>(cfg.getVersion()) < FPGAConfigurationEagleImp::minor_required) {
                stringstream ss;
                ss << "Found FPGA version " << get<0>(cfg.getVersion()) << "." << get<1>(cfg.getVersion()) << "rev" << get<2>(cfg.getVersion())
                   << " which is incompatible with this application. Requiring at least version "
                   << FPGAConfigurationEagleImp::major_required << "." << FPGAConfigurationEagleImp::minor_required << ".";
                StatusFile::addError(ss.str());
                exit(EXIT_FAILURE);
            }
        }

        // when using more than one FPGA is implemented: assert equality of all FPGA configurations
//        if(!all_of(begin(configs) + 1, end(configs), [ &configs ](const FPGAConfigurationEagleImp& other) -> bool { return configs[0] == other; })) {
//            StatusFile::addError("This program requires all used FPGAs to have consistent configurations.");
//            exit(EXIT_FAILURE);
//        }

    } // END if (usefpga)

    if(debug) {
        if (!usefpga)
            cout << "Using no FPGAs";
        else {
            cout << "Using " << hs.getFPGAs().size() << " FPGAs (indices:";
            for(const hybridsys::FPGA& f: hs.getFPGAs())
                cout << " " << f.getIndex();
            cout << ")";
        }
        cout << " and ";
        if (!usegpu)
            cout << "no GPUS.";
        else {
            cout << hs.getGPUs().size() << " GPUs (indices:";
            for(const hybridsys::GPU& g: hs.getGPUs())
                cout << " " << g.getIndex();
            cout << ").";
        }
        cout << endl;
        for(unsigned i = 0; i < configs.size(); i++) { // if using no FPGAs, "configs" is empty.
            cout << "FPGA " << hs.getFPGA(i).getIndex() << " Configuration: " << endl;
            FPGAConfigurationEagleImp c = configs[i];
            cout << c;
        }
    }

    return configs;

}
