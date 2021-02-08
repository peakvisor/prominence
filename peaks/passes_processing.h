#pragma once

#include "peak.h"
#include "common.h"
#include "index.h"
#include "logging.hpp"
#include "distance.h"

#include <string>
#include <deque>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <utility>
#include <algorithm>

struct IndexDistance {
    size_t index;
    double distance;
};

double distanceToContour(Geodesic from, Index radius) {
    PlanarIndex planar{Index(from.latitude), Index(from.longitude)};
    Distance edgeMinDistance = kEarthRadius * 4;
    std::array<Geodesic, 2> region;
    region[0] = Geodesic{double(planar.x - radius), double(planar.y - radius)};
    region[1] = Geodesic{planar.x + radius + 1., planar.y + radius + 1.};
    for (int i = 0; i < 2; ++i) { // const coordinate
        for (int j = 0; j < 2; ++j) { // i value
            for (int k = 0; k <= 2 * radius; ++k) {
                Geodesic sideFrom;
                Geodesic sideTo;
                sideFrom.c[i] = region[j].c[i];
                sideTo.c[i] = region[j].c[i];
                sideFrom.c[1 - i] = region[0].c[1 - i] + k;
                sideTo.c[1 - i] = sideFrom.c[1 - i] + 1;
                sideFrom.normalize();
                sideTo.normalize();
                if (sideFrom.longitude == 180 && sideTo.longitude != 180) {
                    sideFrom.longitude = -180;
                }
                if (sideFrom.latitude == sideTo.latitude && std::abs(sideFrom.latitude) == 90) {
                    continue;
                }
                
                if (std::abs(sideFrom.latitude) == 90) {
                    sideFrom.longitude = sideTo.longitude;
                }
                if (std::abs(sideTo.latitude) == 90) {
                    sideTo.longitude = sideFrom.longitude;
                }
                
                if (sideTo.latitude <= sideFrom.latitude && sideTo.longitude <= sideFrom.longitude) {
                    std::swap(sideTo, sideFrom);
                }
                edgeMinDistance = std::min(DistanceTools::minMaxDistanceToQuadrant(from, sideFrom, sideTo).min, edgeMinDistance);

            }
        }
    }
    return edgeMinDistance;
}

struct Pass {
    uint64_t id;
    std::string name;
    double latitude;
    double longitude;
    Altitude altitude;
    
    void print() {
        std::cout << id << " " << name << " " << latitude << " " << longitude << " " << altitude << std::endl;
    }
};

std::deque<Pass> loadPasses(const std::string &inputDatabase) {
    std::deque<Pass> passes;
    std::string nextLine, nextComponent;
    std::deque<std::string> components;
    std::istringstream scanner;
    
    std::ifstream input{inputDatabase};
    while (std::getline(input, nextLine)) {
        scanner.str(nextLine);
        scanner.clear();
        components.clear();
        while (std::getline(scanner, nextComponent, '\t')) {
            components.push_back(nextComponent);
        }
        if (components.size() < 5) {
            continue;
        }
        
        Pass pass;
        
        pass.id = std::stoull(components[0]);
        pass.name = components[1];
        pass.latitude = std::stod(components[2]);
        pass.longitude = std::stod(components[3]);
        pass.altitude = std::stoi(components[4]);
        
        passes.push_back(pass);
    }
    if (passes.empty()) {
        throw std::logic_error("loaded no passes from " + inputDatabase);
    }
    return passes;
}

struct PeaksAndPasses {
    std::deque<Pass> passes;
    std::deque<Peak> peaks;
};

std::deque<IndexDistance> peaksWithCloseKeycols(const std::deque<Peak> &peaks, const std::map<PlanarIndex, std::deque<size_t>> &sortedOutPeaks, const Geodesic &from, const Altitude &fromAltitude) {
    const double kDistCutoff = 0.3;
    const double kDistEps = 0.01;

    double distanceToTheGrid = distanceToContour(from, 0);
    Index r = distanceToTheGrid > kDistCutoff + kDistEps ? 1 : 0;
    
    std::deque<IndexDistance> result;
    PlanarIndex planar{Index(from.latitude), Index(from.longitude)};
    
    for (Index i = -r; i <= r; ++i) {
        for (Index j = -r; j <= r; ++j) {
            PlanarIndex p = planar + PlanarIndex{i, j};
            auto it = sortedOutPeaks.find(p);
            if (it == sortedOutPeaks.end()) {
                continue;
            }
            for (const size_t &peakIndex : it->second) {
                const auto &peak = peaks[peakIndex];
                assert(peak.keyColLatitude > -1000 && peak.keyColLongitude > -1000);
                auto dist = DistanceTools::distanceOnEarth(from, {peak.keyColLatitude, peak.keyColLongitude});
                if (dist < kDistCutoff) {
                    result.push_back({peakIndex, dist});
                }
            }
        }
    }
    return result;
}

struct PeakWithPasses {
    size_t peakIndex;
    std::deque<IndexDistance> passes;
};

void addPasses(std::deque<Peak> &peaks, const std::string &passesDatabase) {
    const auto passes = loadPasses(passesDatabase);

    std::map<PlanarIndex, std::deque<size_t>> sortedOutPasses;
    std::map<PlanarIndex, std::deque<size_t>> sortedOutPeaks;

    for (size_t i = 0; i < peaks.size(); ++i) {
        const auto &peak = peaks[i];
        if (peak.keyColLatitude == kNoDegrees || peak.keyColLongitude == kNoDegrees) {
            continue;
        }
        assert(peak.keyColLatitude > -1000 && peak.keyColLongitude > -1000);
        PlanarIndex planar{Index(peak.keyColLatitude), Index(peak.keyColLongitude)};
        auto it = sortedOutPeaks.find(planar);
        if (it == sortedOutPeaks.end()) {
            it = sortedOutPeaks.insert({planar, {}}).first;
        }
        it->second.push_back(i);
    }

    for (size_t i = 0; i < passes.size(); ++i) {
        const auto &pass = passes[i];
        PlanarIndex planar{Index(pass.latitude), Index(pass.longitude)};
        auto it = sortedOutPasses.find(planar);
        if (it == sortedOutPasses.end()) {
            it = sortedOutPasses.insert({planar, {}}).first;
        }
        it->second.push_back(i);
    }

    std::map<size_t, std::deque<IndexDistance>> passesWithCloseKeyCols;
    std::set<size_t> passesWithKeycols;

    for (size_t i = 0; i < passes.size(); ++i) {
        const auto &pass = passes[i];
        auto closePeaks = peaksWithCloseKeycols(peaks, sortedOutPeaks, {pass.latitude, pass.longitude}, pass.altitude);
        if (!closePeaks.empty()) {
            passesWithKeycols.insert(i);
        }
        for (const auto &peakIndexDist : closePeaks) {
            auto it = passesWithCloseKeyCols.find(peakIndexDist.index);
            if (it == passesWithCloseKeyCols.end()) {
                it = passesWithCloseKeyCols.insert({i, {}}).first;
            }
            it->second.push_back(peakIndexDist);
        }
    }
    std::map<size_t, std::deque<IndexDistance>> peaksWithPassesMap;
    std::set<size_t> whorePasses;

    for (const auto &passPeaks : passesWithCloseKeyCols) {
        if (passPeaks.second.size() > 1) {
            whorePasses.insert(passPeaks.first);
        }

        for (const auto &peakDistance : passPeaks.second) {
            auto it = peaksWithPassesMap.find(peakDistance.index);
            if (it == peaksWithPassesMap.end()) {
                it = peaksWithPassesMap.insert({peakDistance.index, {}}).first;
            }
            it->second.push_back({passPeaks.first, peakDistance.distance});
        }
    }

    std::deque<PeakWithPasses> peaksWithPasses;
    for (const auto &peakWithPasses : peaksWithPassesMap) {
        peaksWithPasses.push_back({peakWithPasses.first, peakWithPasses.second});
    }

    std::sort(peaksWithPasses.begin(), peaksWithPasses.end(), [&](const PeakWithPasses &first, const PeakWithPasses &second) {
        return peaks[first.peakIndex].prominence > peaks[second.peakIndex].prominence;
    });

    flog() << "total peaksWithPasses: " << peaksWithPasses.size() << std::endl;

    for (const auto &peakWithPasses : peaksWithPasses) {
        if (peakWithPasses.passes.size() != 1) {
            continue;
        }
        const auto passIndex = peakWithPasses.passes.front().index;
        if (whorePasses.count(passIndex)) {
            continue;
        }
        const auto &pass = passes[passIndex];
        auto &peak = peaks[peakWithPasses.peakIndex];
        peak.passId = pass.id;
    }
}

