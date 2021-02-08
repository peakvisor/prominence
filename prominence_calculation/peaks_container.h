#pragma once

#include "common.h"
#include "index.h"
#include "seed.h"
#include "distance.h"
#include "logging.hpp"
#include "dem_grid.h"
#include "peak.h"
#include "../Lib/converter.h"

#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include <cmath>
#include <set>


class PeaksContainer {
 private:
    struct ProminenceDebugInfo {
        Altitude prominence;
        std::string peakName;
        std::string parentName;
    };

    const Converter &converter;
    std::ofstream output;
    
    std::deque<Peak> peaks;
    std::vector<ProminenceDebugInfo> prominenceMapping;
    std::vector<CollectedSeedInfo> seedsInfo;
    std::set<PeakId> mustComputeLeft;
    uint64_t computedProminences = 0;
    
 public:
    PeaksContainer(const Converter &converter, std::deque<Peak> &&peaks,
            const std::string &outputPeaksFile)
            : converter(converter)
            , output(outputPeaksFile)
            , peaks(peaks) { // xcode indentation bug
    }

    template <typename View>
    void init(const DEMGrid<View> &region) {
        output << Peak::header() << std::endl;
        syncWithView(region);
        struct PeakComp { // TODO: remove
            bool operator()(const Peak &f, const Peak &s) const {
                return f.altitude == s.altitude ? f.id < s.id : f.altitude > s.altitude;
            }
        };

        std::sort(peaks.begin(), peaks.end(), PeakComp{});
        
        for (auto &peak : peaks) {
            assert(region.geodesicIsInside(peak.geo()));
            auto planar = region.geodesicToPlanar(peak.geo());
            peak.islandParent = kNoId;
            peak.prominence = kNoAltitude;
            seeds.push_back({{planar, peak.altitude}, kNoIsland});
        }
                
        prepareSeeds();
        for (int i = 1; i < static_cast<int>(seeds.size()); ++i) {
            assert(seeds[i].altitude() <= seeds[i-1].altitude());
            assert(seedsInfo[i].originalIndex == i);
        }
    }
    
    template <typename View>
    void syncWithView(const DEMGrid<View> &region) {
        for (auto &peak : peaks) {
            auto planar = region.geodesicToPlanar({peak.latitude, peak.longitude});
            Altitude fieldAltitude = region.get(planar);
            peak.altitude = fieldAltitude;
            if (!region.isPeak(planar, fieldAltitude)) {
                flog() << "not local maximum: " << peak.id << std::endl;
            }
            if (fieldAltitude > peak.altitude) {
                flog() << "lifting: " << peak.id << " from: " << peak.altitude << " to: " << fieldAltitude << std::endl;
                peak.altitude = fieldAltitude;
            }
            if (fieldAltitude <= kOceanFloor) {
                flog() << "field is too low for: " << peak.id << " a: " << fieldAltitude << std::endl; // BTN
                peak.altitude = 0;
            }
        }
    }
    
    void prepareSeeds() {
        auto renumberSeeds = [&]() {
            int i = 0;
            for (auto &seed : seeds) { seed.islandOwner = i++; }
        };
        
        renumberSeeds();
        
        seedsInfo.clear();

        auto seedsSize = static_cast<IslandIndex>(seeds.size());

        for (auto &seed : seeds) {
            assert(seed.islandOwner < seedsSize);
            seedsInfo.emplace_back(seed.islandOwner);
            if (seed.altitude() <= kOceanFloor) {
                const auto &peak = getPeak(IslandIndex(seedsInfo.size() - 1));
                flog() << peak.name << " ("
                    << peak.id
                << ") at " << Geodesic{peak.latitude, peak.longitude} << " ~ " << seed.planar() << " orig alt: " << peak.altitude << " is low\n";
                seedsInfo.back().reliableProminence = false;
            }
        }
        renumberSeeds();
    }
        
    std::string getPeakName(IslandIndex peakIndex) const {
        return peaks[seedsInfo[peakIndex].originalIndex].name;
    }
    
    Peak& getPeak(IslandIndex peak) {
        return peaks[seedsInfo[peak].originalIndex];
    }
    
    std::string getIslandOwnerName(IslandIndex peak) {
        return peaks[seedsInfo[peak].originalIndex].name;
    }
    
    bool hasReliableProminence(IslandIndex island) {
        if (island == kNoIsland) {
            return false;
        }
        return seedsInfo[island].reliableProminence;
    }
    
    void postProcess() {
        if (failed) {
            return;
        }

        auto seedsSize = static_cast<IslandIndex>(seeds.size());
        for (IslandIndex i = 0; i < seedsSize; ++i) {
            if (seeds[i].islandOwner == i && seedsInfo[i].reliableProminence) {
                auto &peak = getPeak(i);
                peak.prominence = peak.altitude;
                peak.islandParent = peak.id;
                prominenceMapping.emplace_back(ProminenceDebugInfo{peak.prominence, peak.name, ""});
                if (seedsInfo[i].outputed) {
                    flog() << "Already outputed unowned peaks: " << peak.id << std::endl;
                }
            }
        }
        outputRest();
        seedsInfo.clear();
        seeds.clear();
    }
    
    void switchPeakForExtra(const IslandIndex &island, Peak &extraPeak, const SpatialIndex &keycol) {
        auto &switchedPeak = getPeak(island);
        assert(switchedPeak.islandParent == kNoId);
        assert(extraPeak.altitude >= switchedPeak.altitude);
        flog() << "Switching peaks: " << switchedPeak.id << " for: " << extraPeak.id << " island: " << island << std::endl;
        switchedPeak.ilpLatitude = kNoDegrees;
        switchedPeak.ilpLongitude = kNoDegrees;
        extraPeak.islandParent = kNoId;
        Geodesic keycolGeo = converter.planarToGeodesic(keycol.planar);
        setIslandParent(switchedPeak, extraPeak, keycolGeo, keycol.altitude);
        seeds[island].spatial = {converter.geodesicToPlanar({extraPeak.latitude, extraPeak.longitude}), extraPeak.altitude};
        std::swap(switchedPeak, extraPeak);
        assert(seeds[island].spatial.altitude == getPeak(island).altitude);
    }
    
    void addExtraPeaks(std::deque<Peak> &extraPeaks, const SpatialIndex &keycol, IslandIndex parentIsland) {
        parentIsland = getIslandOwner(parentIsland);
        if (!seedsInfo[parentIsland].reliableProminence) {
            return;
        }
        assert(!extraPeaks.empty());
        Altitude highestExtraAltitude = extraPeaks.front().altitude;
        size_t highestIndex = 0;
        for (size_t i = 1; i < extraPeaks.size(); ++i) {
            if (extraPeaks[i].altitude > highestExtraAltitude) {
                highestIndex = i;
                highestExtraAltitude = extraPeaks[i].altitude;
            }
        }
        
        auto &parentPeak = getPeak(parentIsland);
//        if (extraPeaks[highestIndex].id == dbgId) {
//            std::cout << "dbgId is highest\n";
//            std::cout << highestExtraAltitude << " parent: " << parentPeak.toDescriptiveString() << std::endl;
//        }
        if (parentPeak.definitiveAltitudeLess(extraPeaks.at(highestIndex))) {
            auto &highestExtraPeak = extraPeaks[highestIndex];
            auto copyHighest = highestExtraPeak;
            auto copyParent = parentPeak;
            switchPeakForExtra(parentIsland, highestExtraPeak, keycol);
            assert(copyHighest.id == parentPeak.id);
            assert(copyParent.id == highestExtraPeak.id);
        }
        
        unsigned short topPeaks = 0;
        
        for (auto &peak : extraPeaks) {
            if (peak.id == peak.islandParent) {
                assert(++topPeaks < 4);
                setIslandParent(peak, parentPeak, converter.planarToGeodesic(keycol.planar), keycol.altitude);
            }
        }
        
        for (const auto &peak : extraPeaks) {
//            std::cout << "extra: " << peaks.toDescriptiveString() << std::endl;
            output << peak.toStr() << std::endl;
        }
        peaks.insert(peaks.end(), extraPeaks.begin(), extraPeaks.end());
    }
    
    void d__outputProminenceMapping() {
        std::sort(prominenceMapping.begin(), prominenceMapping.end(), [](const ProminenceDebugInfo &one, const ProminenceDebugInfo &other) {
            return one.prominence > other.prominence;
        });
        static const int kPartOfInterest = 10;
        flog() << "Top " << kPartOfInterest << " peaks by prominence:\n";
        for (int i = 0; i < kPartOfInterest && i < static_cast<int>(prominenceMapping.size()); ++i) {
            const auto &info = prominenceMapping[i];
            flog() << "\t" << info.prominence << "\t" << info.peakName << "\t->\t" << info.parentName <<  "\n";
        }
    }
        
    IslandIndex getIslandOwner(IslandIndex island) {
        IslandIndex nextOwner = island, owner, chainLength = 0;
        do {
            ++chainLength;
            owner = nextOwner;
            nextOwner = seeds[owner].islandOwner;
        } while (nextOwner != owner);
        
        int iter = island;
        while (--chainLength > 0) {
            auto &seed = seeds[iter];
            iter = seed.islandOwner;
            seed.islandOwner = owner;
        }
        return owner;
    }
    
    IslandIndex subductIslands(IslandIndex one, IslandIndex another, const SpatialIndex &col) {
        one = getIslandOwner(one);
        another = getIslandOwner(another);
        if (one == another) {
            return one;
        }
        
        bool anotherIsDominant = getPeak(one).definitiveAltitudeLess(getPeak(another));
        auto dominant = anotherIsDominant ? another : one;
        auto inferior = anotherIsDominant ? one : another;
        
        auto &inferiorInfo = seedsInfo[inferior];
        auto const &dominantInfo = seedsInfo[dominant];
        auto &inferiorPeak = getPeak(inferior);
        auto const &dominantPeak = getPeak(dominant);
        seeds[inferior].islandOwner = dominant;
        
        if (inferiorPeak.isExtra()) {
            flog() << "Subduction inferior: " << inferior << " id: "  << inferiorPeak.id << "(reliable: " << inferiorInfo.reliableProminence << "|" << dominantInfo.reliableProminence  << ") dominant: " << dominant << " id:" << dominantPeak.id << std::endl;
        }
        Geodesic colGeo = converter.planarToGeodesic(col.planar);
        if (inferiorPeak.altitude != seeds[inferior].altitude()) {
            std::cout << inferiorPeak.toDescriptiveString() << std::endl;
            std::cout << dominantPeak.toDescriptiveString() << std::endl;
            assert(false);
        }
        
        if (inferiorInfo.reliableProminence && dominantInfo.reliableProminence) {
            setIslandParent(inferiorPeak, dominantPeak, colGeo, col.altitude);
            prominenceMapping.emplace_back(ProminenceDebugInfo{inferiorPeak.prominence, inferiorPeak.name, dominantPeak.name});
 
            if (inferiorInfo.outputed) {
                flog() << "err: trying to output " << inferiorPeak.id << " second time\n";
            } else {
                assert(inferiorPeak.islandParent != kNoId);
                outputPeak(inferior);
            }
        } else {
            if (inferiorInfo.reliableProminence) {
                if (inferiorInfo.outputed) {
                    flog() << "2err: trying to output " << inferiorPeak.id << " second time\n";
                } else {
                    setProminence(inferiorPeak, colGeo, col.altitude);
                    outputPeak(inferior);
                }
            }
            if (dominantInfo.reliableProminence) {
                markUnreliableProminence(dominant);
            }
        }
        return dominant;
    }
    
    bool markUnreliableProminence(IslandIndex island) {
        IslandIndex actualOwner = getIslandOwner(island);

        auto &info = seedsInfo[actualOwner];
        if (info.reliableProminence) {
            flog() << "marking " << getPeak(actualOwner).id << std::endl;
            info.reliableProminence = false;
            if (info.outputed) {
                flog() << "err: trying to output " << getPeak(actualOwner).id << " second time while marking unreliable\n";
            } else {
                auto &peak = getPeak(actualOwner);
                peak.prominence = kNoAltitude;
                outputPeak(island);
            }
            return true;
        }
        return false;
    }
    
    void setIsolation(IslandIndex islandOwner, Isolation &&isolation) {
        auto &peak = getPeak(islandOwner);
        Geodesic ilp = converter.planarToGeodesic(isolation.ilp);
        peak.ilpLatitude = ilp.latitude;
        peak.ilpLongitude = ilp.longitude;
    }
    
    void outputRest() {
        auto size = static_cast<IslandIndex>(seedsInfo.size());
        for (IslandIndex i = 0; i < size; ++i) {
            if (!seedsInfo[i].outputed) {
                outputPeak(i);
            }
        }
    }
    
    void markAllLeftUnreliable() {
        for (auto &info : seedsInfo) {
            if (!info.outputed) {
                info.reliableProminence = false;
            }
        }
    }
    
    bool finished() {
        return computedProminences == seeds.size();
    }
    
    void outputPeak(IslandIndex island) {
        assert(!seedsInfo[island].outputed);
        ++computedProminences;
        seedsInfo[island].outputed = true;
        auto &peak = getPeak(island);
        output << peak.toStr() << std::endl;
        if (!output) {
            flog() << "err: broken output\n";
            std::cout << "err: broken output\n";
            assert(false);
        }
        assert(output);
    }
    
    std::deque<Peak> unloadPeaks() {
        std::deque<Peak> result;
        std::swap(peaks, result);
        return result;
    }
    
    std::vector<Seed> seeds;
    bool failed = false;
};

