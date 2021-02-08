#pragma once

#include "index.h"

#include <functional>

struct Isolation {
    double distance;
    PlanarIndex ilp;
    IslandIndex isolationParent = kNoIsland;
    bool reliable;
};

struct CollectedSeedInfo {
    const int originalIndex;
    bool reliableProminence = true;
    bool outputed = false;

    CollectedSeedInfo(int originalIndex) :
        originalIndex(originalIndex) {}
};

struct Seed {
    SpatialIndex spatial;
    mutable IslandIndex islandOwner;
    
    Altitude altitude() const {
        return spatial.altitude;
    }
    
    const PlanarIndex &planar() const {
        return spatial.planar;
    }
    
    inline bool operator<(const Seed &other) const {
        if (spatial == other.spatial) {
            if (islandOwner == kNoIsland && other.islandOwner != kNoIsland) {
                return true;
            }
            if (other.islandOwner == kNoIsland && islandOwner != kNoIsland) {
                return false;
            }
            return islandOwner < other.islandOwner;
        }
        return spatial < other.spatial;
    }
    
//    bool operator >(const Seed &other) const {
//        return altitude == other.altitude ? planar > other.planar : altitude > other.altitude;
//    }
    
//    inline bool operator ==(const Seed &other) const {
//        return planar.prehash == other.planar.prehash; //x == other.x && y == other.y;
//    }
};

struct AdditionalPeak {
    SpatialIndex spatial;
    SpatialIndex keyColLike;
    IslandIndex parentLike = kNoIsland;
    uint64_t area{0};
    uint64_t volume{0};
    PlanarIndex lastSeenPlanar{0, 0};
    
    Altitude prominence() const {
        return spatial.altitude - keyColLike.altitude;
    }
    
    bool isInited() const {
        return parentLike != kNoIsland;
    }
    
    void deInit() {
        parentLike = kNoIsland;
    }
};
