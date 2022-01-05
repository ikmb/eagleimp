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

#include <cstdlib>
#include <string>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>

#include "MapInterpolater.h"

#include "StatusFile.h"

using namespace std;

const string MapInterpolater::MAP_FILE_HEADER1 =
        "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)";
const string MapInterpolater::MAP_FILE_HEADER2 =
        "position COMBINED_rate(cM/Mb) Genetic_Map(cM)";

MapInterpolater::MapInterpolater(const string &geneticMapFile, int pick_chr_) : pick_chr(pick_chr_) {
    cout << "Reading genetic map coordinates of chromosome " << pick_chr << "." << endl;
    chrBpToRateGen[0] = make_pair(0.0, 0.0); // sentinel at beginning

    int bp0 = 0;
    double gen0 = 0;
    int bp;
    double rate, gen;

    // normal case (not chrY)
    if (pick_chr != CHRY) {

        string line;
        ifstream fin;
        fin.open(geneticMapFile);
        if (!fin.good()) {
            StatusFile::addError("Failed reading genetic map.");
            exit(EXIT_FAILURE);
        }
        getline(fin, line);
        if (line != MAP_FILE_HEADER1 && line != MAP_FILE_HEADER2) {
            StatusFile::addError("Wrong format of genetic map.\n"
            "       Expecting header either: " + MAP_FILE_HEADER1 + "\n"
            "                            or: " + MAP_FILE_HEADER2);
            exit(EXIT_FAILURE);
        }
        if (line == MAP_FILE_HEADER1) { // map contains chromosome column
            int chr;
            bool reachedpick = false;
            while (fin >> chr >> bp >> rate >> gen) {
                if (reachedpick) {
                    if (chr == pick_chr) {
                        if (gen < gen0 || bp < bp0) {
                            StatusFile::addError("Genetic map contains out-of-order row.");
                            exit(EXIT_FAILURE);
                        }
                        if (bp <= 0) {
                            StatusFile::addError("Genetic map positions must be positive.");
                            exit(EXIT_FAILURE);
                        }
                        if (bp == bp0) {
                            StatusFile::addError("Genetic map contains duplicate position.");
                            exit(EXIT_FAILURE);
                        }
                        chrBpToRateGen[bp] = make_pair((gen - gen0) / (1e-6 * (bp - bp0)), gen);
                    } else // all of the chromosome was consumed, we can stop here
                        break;
                } else {
                    if (chr == pick_chr)
                        reachedpick = true;
                }
                bp0 = bp;
                gen0 = gen;
            }
        } else { // no chromosome column in map -> we expect everything to be on the right chromosome!
            while (fin >> bp >> rate >> gen) {
                if (gen < gen0 || bp < bp0) {
                    StatusFile::addError("Genetic map contains out-of-order row.");
                    exit(EXIT_FAILURE);
                }
                if (bp <= 0) {
                    StatusFile::addError("Genetic map positions must be positive.");
                    exit(EXIT_FAILURE);
                }
                if (bp == bp0) {
                    StatusFile::addError("Genetic map contains duplicate position.");
                    exit(EXIT_FAILURE);
                }
                chrBpToRateGen[bp] = make_pair((gen - gen0) / (1e-6 * (bp - bp0)), gen);
                bp0 = bp;
                gen0 = gen;
            }
        }

    }
    chrBpToRateGen[bp0 + 1e9] = make_pair(1.0, gen0 + 1e3); // sentinel at end -> this also generates equidistant positions with 1cM/Mbp for special case with chrY (where no phasing is required)
}

// returns interpolated genetic position in cM (centi morgans)
double MapInterpolater::interp(int bp) const {
    map<int, pair<double, double> >::const_iterator ubIter =
            chrBpToRateGen.upper_bound(bp); // map record > bp
    if (ubIter == chrBpToRateGen.end()) {
        // better extrapolate! Then, we don't need a sentinel.
        StatusFile::addError("Marker at " + to_string(bp) + " is larger than any in genetic map.");
        exit(EXIT_FAILURE);
    }
    int ubBp = ubIter->first;
    double ubRate = ubIter->second.first;
    double ubGen = ubIter->second.second;

    return ubGen + 1e-6 * ubRate * (bp - ubBp); // interpolate interval
}
