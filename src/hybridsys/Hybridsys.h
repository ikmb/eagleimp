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

#ifndef HYBRIDSYS_H
#define HYBRIDSYS_H

#include <vector>

#include "FPGA.h"
#include "GPU.h"
#include "Buffer.h"
#include "BufferFactory.h"

namespace hybridsys {

/**
 * @brief The Hybridsys class represents the whole system and is meant to be used as the main entry point for all GPU/FPGA activities.
 */
class Hybridsys
{
public:

    /**
     * @brief Hybridsys Initializes and discovers all necessary APIs.
     *
     * The default constructor performs device discovery, i.e. the FPGA and GPU drivers are queried for supported
     * devices. This may take several seconds. Additionally, all discovered FPGAs are locked so they cannot be
     * used by other programs or instances of Hybridsys.
     */
    Hybridsys(std::vector<unsigned> allowed_fpgas = {}, std::vector<unsigned> allowed_gpus = {});

    /**
     * @brief getFPGAs returns a list of usable FPGAs.
     *
     * Returns a list of all usable FPGAs.
     * @note Although all FPGAs delivered via this method are guaranteed to be actually available, they may not
     *  have been configured with a PCIe tandem stage 2 bitstream, and therefore, may not be usable until properly programmed. See also: FPGA::isConfigured()
     */
    std::vector<FPGA>& getFPGAs() { return fpgas; }
    FPGA& getFPGA(unsigned i) { return fpgas.at(i); }
    std::vector<GPU>& getGPUs() { return gpus; }
    GPU& getGPU(unsigned i) { return gpus.at(i); }

private:
    std::vector<FPGA> fpgas;
    std::vector<GPU> gpus;
};

}
#endif // HYBRIDSYS_H
