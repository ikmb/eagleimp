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

#include <vector>
#include <sstream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <iostream>

extern "C" {
#include <sys/types.h>
#include <dirent.h>
}

#include "BufferFactory.h"

#include "Hybridsys.h"

using namespace std;

namespace hybridsys {

struct pci_ident {
    int vendor;
    int device;
};

struct pci_device {
    int bus;
    int slot;
};

static std::vector<struct pci_ident> pci_ident_fpga = {
    {0x4144, 0xADB3}, // Alpha Data; ADM-PCIE-8K5
};
static std::vector<struct pci_ident> pci_ident_gpu = {
    {0x10DE, 0x15F8}, // NVIDIA; Tesla P100 16G/PCIE
    {0x10DE, 0x1DB6}, // NVIDIA; Tesla V100 32G/PCIE
    {0x10DE, 0x137a}  // Nvidia Quadro something...
    // TODO this list needs to be extended to contain all capable Nvidia GPUs, OR: probably only compare vendor ID?
};

#if defined(USE_CUDA_GPU) || defined(USE_AD_FPGA)

static bool operator==(const struct pci_ident &lhs, const struct pci_ident &rhs) {
    return (lhs.vendor == rhs.vendor) && (lhs.device == rhs.device);
}

static int readHex(const std::string& file) {
    int result;
    std::ifstream f(file);
    f >> std::hex >> result;
    f.close();
    return result;
}

static std::vector<struct pci_device> findDevices(const std::vector<struct pci_ident>& idents) {
    std::vector<struct pci_device> devices;

    DIR* dir = opendir("/sys/bus/pci/devices");
    struct dirent *entry = nullptr;

    while((entry = readdir(dir))) {
        if(entry->d_type == DT_LNK) {
            struct pci_ident this_device;
            std::stringstream ss;
            ss << "/sys/bus/pci/devices/";
            ss << entry->d_name;
            ss << "/vendor";

            this_device.vendor = readHex(ss.str());
            ss = std::stringstream();
            ss << "/sys/bus/pci/devices/";
            ss << entry->d_name;
            ss << "/device";
            this_device.device = readHex(ss.str());

            if(std::find(std::begin(idents), std::end(idents), this_device) != std::end(idents)) {
                int bus = 0, slot = 0;
                std::string token;
                std::istringstream name(entry->d_name);

                std::getline(name, token, ':'); // extract domain and skip it
                name >> std::hex >> bus;
                std::getline(name, token, ':'); // extract bus remainder and skip it
                name >> std::hex >> slot;


                devices.push_back( {bus, slot} );
            }
        }
    }
    return devices;
}
#endif

Hybridsys::Hybridsys(std::vector<unsigned> allowed_fpgas, std::vector<unsigned> allowed_gpus)
    : fpgas(),
      gpus()
{
    std::vector<struct pci_device> devices;

    // Alpha Data API does not need to be initialized
#ifdef USE_AD_FPGA
    devices = findDevices(pci_ident_fpga);


    for(const pci_device& d: devices) {
        if(find(begin(allowed_fpgas), end(allowed_fpgas), FPGA::findIndex(d.bus, d.slot)) != end(allowed_fpgas))
            fpgas.emplace_back(d.bus, d.slot);
    }
#else
    (void) allowed_fpgas;
#endif

#ifdef USE_CUDA_GPU
    // Initialize NVML
    nvmlReturn_t ret = nvmlInit();
    if (ret != NVML_SUCCESS)
        throw std::runtime_error(nvmlErrorString(ret));
    devices.clear();
    devices = findDevices(pci_ident_gpu);

    // Filter GPUs
    for(const pci_device& d: devices) {
    	gpus.emplace_back(d.bus, d.slot);
    	if(find(begin(allowed_gpus), end(allowed_gpus), gpus.back().getIndex()) == end(allowed_gpus))
    		// element is not found in the allowed list -> erase
    		gpus.pop_back();
    }
#else
    (void) allowed_gpus;
#endif
}


}
