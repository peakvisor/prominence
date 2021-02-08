#pragma once

#include "common.h"

#include <deque>
#include <set>
#include <functional>
#include <map>
#include <fstream>

struct Peak {
    PeakId id{};
    std::string name{};
    Degrees latitude{kNoDegrees};
    Degrees longitude{kNoDegrees};
    Altitude altitude{kNoAltitude};
    Altitude prominence{kNoAltitude};
    Altitude proportionalProminence{kNoAltitude};
    PeakId islandParent{};
    Degrees keyColLatitude{kNoDegrees};
    Degrees keyColLongitude{kNoDegrees};
    Degrees ilpLatitude{kNoDegrees};
    Degrees ilpLongitude{kNoDegrees};
    PeakId nhn{};
    PeakId passId{};

    static std::string header(); // TODO: generalize

    bool isExtra() const { return name.empty(); }
    Geodesic geo() const { return {latitude, longitude}; }
    Geodesic keyColGeo() const { return {keyColLatitude, keyColLongitude}; }
    Geodesic ilp() const { return {ilpLatitude, ilpLongitude}; }
    std::string toDescriptiveString() const;
    std::string toStr() const; // TODO: generalize
    bool isOrphan() const { return islandParent == kNoId; }

    [[nodiscard]] Altitude keyColAltitude() const {
        return prominence == kNoAltitude ? kOceanFloor : altitude - prominence;
    }

    [[nodiscard]] bool definitiveKeyColAltitudeLess(const Peak &other) const {
        auto keyColAlt = keyColAltitude();
        auto otherKeyColAlt = other.keyColAltitude();
        if (keyColAlt == otherKeyColAlt) {
            if (keyColGeo() == other.keyColGeo()) {
                return id < other.id;
            }
            return keyColGeo() < other.keyColGeo();
        }
        return keyColAlt < otherKeyColAlt;
    }

    [[nodiscard]] bool definitiveAltitudeLess(const Peak &other, bool lowerExtras = true) const {
        if (altitude == other.altitude) {
            if (lowerExtras) {
                bool extra = isExtra();
                bool oextra = other.isExtra();
                if (extra != oextra) {
                    return extra;
                }
            }
            if (geo() == other.geo()) {
                return id < other.id;
            }
            return geo() < other.geo();
        } else {
            return altitude < other.altitude;
        }
    }
};

std::deque<Peak> loadPeaks(const Path &inputDatabase,
        const std::function<bool(const Peak&)> &filter = [](const Peak &) { return true;});
std::deque<Peak> loadPeaks(const Path &inputDatabase, Geodesic from, Geodesic to);
void savePeaks(std::deque<Peak> &peaks, const Path &outputPeaksPath);

void correctIsolation(Peak &peak, const Peak &nhn);
PeakId extraPeakId(Geodesic degrees);
void setProminence(Peak &peak, const Geodesic &keycolGeo, const Altitude &keycolAltitude);
void setIslandParent(Peak &peak, const Peak &parentPeak, const Geodesic &keycolGeo, const Altitude &keycolAltitude);
std::map<PeakId, size_t> getIdToIndexMap(const std::deque<Peak> &peaks, bool resolveDuplicateIds = false);
void sortById(std::deque<Peak> &peaks);
void checkSortedById(const std::deque<Peak> &peaks);

enum DualState {
    kDualStateBoth = 0,
    kDualStateFirst,
    kDualStateSecond
};

using DualIterationCallback = std::function<void(size_t, size_t, DualState)>;
void dualIteration(std::deque<Peak> &first, std::deque<Peak> &second,
    const DualIterationCallback &callback, bool assumeSorted = false);

using UpdateCallback = std::function<void(Peak&, const Peak&)>;
void updateSortedPeaks(std::deque<Peak> &peaks, std::deque<Peak> &updates,
    const UpdateCallback &updateCallback);
bool peakFilesAreTheSame(const Path &lhs, const Path &rhs, bool printDiffPeaks = true);