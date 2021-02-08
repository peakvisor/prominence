#pragma once

#include "common.h"
#include "index.h"
#include "dem_grid.h"

struct PeaksCollector {
    const Altitude cutoff;
    
    PeaksCollector(const Altitude &cutoff) : cutoff(cutoff) {
        deinit();
    }
    
    template <typename View>
    void visit(const DEMGrid<View> &region, const SpatialIndex &spatial) {
        if (spatial.altitude > highest.altitude) {
            highest = spatial;
        }
        if (spatial.altitude < minAltitude() || !region.isPeak(spatial.planar, spatial.altitude)) {
            return;
        }
        
        addPeak(region, spatial);
    }
    
    template <typename View>
    void addPeak(const DEMGrid<View> &region, const SpatialIndex &spatial) {
        peaks.emplace_back();
        auto &peak = peaks.back();
        Geodesic geo = region.planarToGeodesic(spatial.planar);
        peak.altitude = spatial.altitude;
        peak.latitude = geo.latitude;
        peak.longitude = geo.longitude;
        peak.id = extraPeakId(geo);
    }
    
    void deinit() {
        peaks.clear();
        inited = false;
        parent = kNoIsland;
//        parentAltitude = kNoAltitude;
    }
    
    void init(const SpatialIndex &keycol, const SpatialIndex &first, const IslandIndex &parentIsland) {
        assert(peaks.empty());
        assert(!inited);
        assert(parentIsland != kNoIsland);
        inited = true;
        highest = first;
        parent = parentIsland;
        this->keycol = keycol;
    }
    
    const SpatialIndex &getHighest() const {
        return highest;
    }
    
    Altitude minAltitude() const  {
        return keycol.altitude + cutoff;
    }
    
    bool inited = false;
    size_t highestPeakIndex = 0;
    SpatialIndex keycol;
    std::deque<Peak> peaks;
    IslandIndex parent = kNoIsland;
    SpatialIndex highest;
//    Altitude parentAltitude = kNoAltitude;
};
