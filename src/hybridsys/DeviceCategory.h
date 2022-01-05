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

#ifndef DEVICECATEGORY_H
#define DEVICECATEGORY_H

#include <system_error>

namespace hybridsys {

/**
 * @brief An error category to encapsulate runtime errors from the native Alpha Data FPGA API
 */
class FPGACategory : public std::error_category
{
public:
    FPGACategory() : std::error_category() {}

    // error_category interface
public:
    const char *name() const noexcept override ;
    std::string message(int code) const override;
};

/**
 * @brief An error category to encapsulate runtime errors from the native CUDA API
 */
class GPUCategory : public std::error_category
{
public:
    GPUCategory() : std::error_category() {}

    // error_category interface
public:
    const char *name() const noexcept override;
    std::string message(int code) const override;
};

}
#endif // DEVICECATEGORY_H
