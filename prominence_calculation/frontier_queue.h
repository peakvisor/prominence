#pragma once
#include "index.h"
#include "seed.h"

#include <functional>
#include <deque>
#include <array>
#include <cassert>
#include <fstream>
#include <vector>

class FrontierHeap {
 public:
    void push(const SpatialIndex &point) {
        frontier.push_back(point);
        std::push_heap(frontier.begin(), frontier.end());
    }
    
    void push(PlanarIndex planar, Altitude altitude) {
        push(SpatialIndex{planar, altitude});
    }
    
    SpatialIndex pop() {
        assert(!frontier.empty());
        std::pop_heap(frontier.begin(), frontier.end());
        auto next = frontier.back();
        frontier.pop_back();
        return next;
    }
    
    bool empty() const {
        return frontier.empty();
    }
    
    Altitude nextAltitude() const {
        return frontier.empty() ? kOceanFloor : frontier.front().altitude;
    }
    
    size_t size() {
        return frontier.size();
    }
    
 private:
    std::vector<SpatialIndex> frontier;
};

class FrontierQueue {
 public:
    FrontierQueue(const std::vector<Seed> &seeds, const Altitude &maxAltitudeToIgnore)
    : maxAltitudeToIgnore(maxAltitudeToIgnore), seeds(seeds), waterLevel(seeds.front().spatial.altitude) {}
    
    bool push(PlanarIndex planar, Altitude altitude) {
        if (altitude <= maxAltitudeToIgnore) {
            return false;
        }
        
        auto altLevel = altitudeLevel(altitude);
        
        switch (altLevel) {
            case kBelow:
                belowFrontier.push({planar, altitude});
                break;
            case kAbove:
                aboveFrontier.push({planar, altitude});
                break;
            default:
                mainFrontier[levelIndex(altLevel)].push_back(planar);
                break;
        }
        return true;
    }
    
    Seed pop() {
        assert(waterLevel > kOceanFloor);
        if (!aboveFrontier.empty()) {
            auto point = aboveFrontier.pop();
            return {point, kNoIsland};
        }
        
        liftWater();
        Seed seed{};
        
        if (nextSeedToVisit < static_cast<int>(seeds.size()) && seeds[nextSeedToVisit].spatial.altitude == waterLevel) {
            seed = seeds[nextSeedToVisit];
            ++nextSeedToVisit;
        } else {
            assert(!mainFrontier[topLevel].empty());
            auto &topLevelFrontier = mainFrontier[topLevel];
            seed.spatial = {topLevelFrontier.back(), waterLevel};
            seed.islandOwner = kNoIsland;
            topLevelFrontier.pop_back();
        }
//        liftWater();
        
        return seed;
    }
    
    Altitude getWaterLevel() const {
        return waterLevel;
    }
    
    int lastPoppedSeed() const {
        return nextSeedToVisit - 1;
    }
    
    void d__printState() {
        std::cout << "WaterLevel: " << waterLevel << " nextSeedAltitude: " << seeds[nextSeedToVisit].altitude()
        << " aboveHeapSize: " << aboveFrontier.size() << " belowHeapSize: " << belowFrontier.size() << std::endl;
        for (int i = 0; i < kLevelsToStore; ++i) {
            std::cout << " on frontier " << i << ": " << mainFrontier[(topLevel + i) % kLevelsToStore].size() << "\n";
        }
    }
    
    void liftWater() {
        if (!aboveFrontier.empty()) {
            return;
        }

        while (mainFrontier[topLevel].empty()
            && (nextSeedToVisit >= static_cast<int>(seeds.size()) || seeds[nextSeedToVisit].altitude() < waterLevel))
        {
            --waterLevel;
            if (waterLevel < kOceanFloor) {
                return;
            }
            topLevel = (topLevel + 1) % kLevelsToStore;
            Altitude bottomLevel = waterLevel - kLevelsToStore + 1;
            if (bottomLevel > kOceanFloor) {
                while (belowFrontier.nextAltitude() == bottomLevel) {
                    mainFrontier[bottomLevelIndex()].push_back(belowFrontier.pop().planar);
                }
            }
        }
    }
    
private:
    static const int kLevelsToStore = 100;
    static const int kAbove = -1;
    static const int kBelow = -2;
    const Altitude maxAltitudeToIgnore;
    const std::vector<Seed> &seeds;
    
    int altitudeLevel(Altitude alt) {
        auto diff = waterLevel - alt;
        if (diff < 0) {
            return kAbove;
        }
        if (diff >= kLevelsToStore) {
            return kBelow;
        } else {
            return diff;
        }
    }
    
    int levelIndex(int altLevel) const {
        return (topLevel + altLevel) % kLevelsToStore;
    }
    
    int bottomLevelIndex() const {
        return (topLevel + (kLevelsToStore - 1)) % kLevelsToStore;
    }
    
    Altitude waterLevel;
    
    FrontierHeap belowFrontier;
    FrontierHeap aboveFrontier;
    std::array<std::deque<PlanarIndex>, kLevelsToStore> mainFrontier;
    
    int topLevel = 0;
    int nextSeedToVisit = 0;
};

class SeedHeap {
 public:
    SeedHeap(const std::vector<Seed> &seeds, const Altitude &maxAltitudeToIgnore) : seeds(seeds), maxAltitudeToIgnore(maxAltitudeToIgnore), nextSeedToVisit(0) {
        waterLevel = nextAltitude();
    }
    
    bool push(const Seed &seed) {
        if (seed.altitude() <= maxAltitudeToIgnore) {
            return false;
        }
        heap.push(seed.spatial);
        return true;
    }
    
    bool push(PlanarIndex planar, Altitude altitude, IslandIndex islandIndex = kNoIsland) {
        return push(Seed{{planar, altitude}, islandIndex});
    }
    
    Seed pop() {
        liftWater();
        if (nextSeedToVisit < seeds.size() && seeds[nextSeedToVisit].altitude() == waterLevel) {
            return seeds[nextSeedToVisit++];
        }
        return {heap.pop(), kNoIsland};
    }
    
    void liftWater() {
        if (nextAltitude() < waterLevel) {
            waterLevel = nextAltitude();
        }
    }
    
    Altitude nextAltitude() const {
        Altitude nextSeedAltitude = nextSeedToVisit < seeds.size() ? seeds[nextSeedToVisit].altitude() : maxAltitudeToIgnore;
        Altitude nextHeapAltitude = heap.empty() ? maxAltitudeToIgnore : heap.nextAltitude();
        return std::max(nextSeedAltitude, nextHeapAltitude);
    }

    Altitude getWaterLevel() const {
        return waterLevel;
    }
    
 private:
//    struct RCmp {
//        bool operator()(const Seed &f, const Seed &s) const {
//            return s < f;
//        };
//    };
    
    const std::vector<Seed> &seeds;
    const Altitude maxAltitudeToIgnore;
    FrontierHeap heap;
    size_t nextSeedToVisit = 0;
    Altitude waterLevel;
};

