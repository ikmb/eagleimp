/*
 *    Copyright (C) 2018-2021 by Lars Wienbrandt,
 *    Institute of Clinical Molecular Biology, Kiel University
 *    and 2015-2016 Po-Ru Loh, Harvard University
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

#ifndef MAPINTERPOLATER_HPP
#define MAPINTERPOLATER_HPP

#include <string>
#include <map>
#include <utility>

using namespace std;

static const int CHRX = 23; // mapping chromosome X to number 23
static const int CHRY = 24; // mapping chromosome Y to number 24

class MapInterpolater {

public:
    MapInterpolater(){}
    // input file format either: chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)
    // (Oxford map format preceded by chr column),
    // then it just picks the chromosome indicated with chr
    // OR position COMBINED_rate(cM/Mb) Genetic_Map(cM)
    // requiring a single chromosome (the target chromosome) in the map
    MapInterpolater(const std::string &geneticMapFile, int chr);
    // returns interpolated genetic position in Morgans
    double interp(int bp) const;

private:
    map<int, std::pair<double, double> > chrBpToRateGen;
    static const string MAP_FILE_HEADER1;
    static const string MAP_FILE_HEADER2;
    int pick_chr = 0;
};

#endif
