#ifndef PROMINENCE_CALCULATOR_PEAKS_ON_DEM_PROCESSING_H
#define PROMINENCE_CALCULATOR_PEAKS_ON_DEM_PROCESSING_H

#include "index.h"
#include "dem_filename_formatter.h"
#include "peak.h"
#include "logging.hpp"
#include "dem_grid.h"

#include <random>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <map>
#include <algorithm>
#include <set>

Altitude spikeSuspicion(const std::string &demsDir, const Peak &peak, Index radius, int diffCutoff) {
    PlanarIndex demPlanar{static_cast<Index>(std::floor(peak.latitude)), static_cast<Index>(std::floor(peak.longitude))};
    TiledDEMGrid region(demsDir, demPlanar, demPlanar, {kBigTiledFileSizeX, kBigTiledFileSizeY}, 3600);
    if (region.filesCount == 0) {
        return -1;
    }

    PlanarIndex peakPlanar = region.geodesicToPlanar(peak.geo());
    assert(region.planarIsInside(peakPlanar));

    Rectangle searchRectangle{peakPlanar - kDiagonalShift * radius, peakPlanar + kDiagonalShift * radius};
    Altitude diff = 0;
    for (Index la = searchRectangle.origin.latitude - 1; la < searchRectangle.end.latitude; ++la) {
        for (Index lo = searchRectangle.end.longitude - 1; lo < searchRectangle.end.longitude; ++lo) {
            if (!region.planarIsInside({la, lo})) {
                continue;
            }
            Altitude alt = region.get({la, lo});
            for (const auto &neighbour : {PlanarIndex{la + 1, lo}, PlanarIndex{la, lo + 1}}) {
                if (!region.planarIsInside(neighbour)) {
                    continue;
                }
                Altitude newDiff = std::abs(alt - region.get(neighbour));
                if (newDiff > diff) {
                    diff = newDiff;
                }
            }
        }
    }
    return diff;
}

void findExtraUltrasCandidates(const std::string &postprocessedDatabase, const std::string &referenceDatabase, const std::string &demsDir) {
    auto peaks = loadPeaks(postprocessedDatabase, [](const Peak &peak) {
        return peak.isExtra() && peak.prominence >= 500;
    });
    auto refpeaks = loadPeaks(referenceDatabase, [](const auto& peak) {
        return peak.isExtra() && peak.prominence >= 1500;
    });
    auto idToIndex = getIdToIndexMap(refpeaks);

    std::vector<std::pair<Altitude, Index>> spikenessIndex;

    for (Index i = 0; i < static_cast<Index>(peaks.size()); ++i) {
        const auto &peak = peaks.at(i);
        if (peak.prominence >= 1500 || idToIndex.find(peak.id) != idToIndex.end()) {
            spikenessIndex.emplace_back(spikeSuspicion(demsDir, peak, 10, 500), i);
        }
    }
    std::sort(spikenessIndex.begin(), spikenessIndex.end());
    for (const auto &p : spikenessIndex) {
        std::cout << "spikeness: " <<  p.first << "\t";
        std::cout << peaks.at(p.second).toDescriptiveString() << std::endl;
    }
}

void applyPeaksOnDemProcessing(int argc, char *argv[]) {
    std::string mode = argc > 1 ? std::string(argv[1]) : "in_code"s;
    if (mode == "find_extra_ultras") {
        findExtraUltrasCandidates(argv[2], argv[3], argv[4]);
    }
}

#endif //PROMINENCE_CALCULATOR_PEAKS_ON_DEM_PROCESSING_H
