#pragma once
#include "logging.hpp"
#include "peak.h"
#include "index.h"

#include <deque>
#include <vector>
#include <map>
#include <algorithm>

struct PeakTree {
    struct Node {
        size_t index = kNoIndexSizeT;
        size_t parent = kNoIndexSizeT;
        std::deque<size_t> islandKids;
        
        size_t kidIndexInParent(size_t kidIndex) {
            for (size_t i = 0; i < islandKids.size(); ++i) {
                if (kidIndex == islandKids[i]) {
                    return i;
                }
            }
            return kNoIndexSizeT;
        }
        
        bool isOrphan() const {
            return parent == kNoIndexSizeT;
        }
        
    };
    
    struct KeyCol {
        size_t parent;
        std::deque<size_t> kids;
    };
    
    std::deque<Peak> &peaks;
    std::deque<Node> nodes;
    std::map<PeakId, size_t> idToIndex;
    std::map<PeakId, KeyCol> keyCols;
    
    explicit PeakTree(std::deque<Peak> &peaks, bool resolveDuplicateIds = false, bool fillKeycols = false) : peaks(peaks) {
        idToIndex = getIdToIndexMap(peaks, resolveDuplicateIds);
        assert(peaks.size() == idToIndex.size());
        buildTree();
        if (fillKeycols) {
            resolveParentageInKeycols();
            fillKeyCols();
        }
    }
    
    void buildTree() {
        nodes.resize(peaks.size());
        for (size_t i = 0; i < peaks.size(); ++i) {
            const auto &peak = peaks[i];
            nodes[i].index = i;
            if (peak.islandParent == peak.id) {
                nodes[i].parent = i;
            } else if (peak.islandParent != kNoId) {
//                std::cout << peaks.id << endl;
                auto it = idToIndex.find(peak.islandParent);
                if (it != idToIndex.end()) {
                    nodes.at(i).parent = it->second;
                } else {
                    std::cout << peaks[i].toDescriptiveString() << " is bogus parent wise\n";
                    assert(false);
                }
                if (peak.islandParent != peak.id) {
                    nodes[nodes[i].parent].islandKids.push_back(i);
                }
            }
        }
        for (size_t i = 0; i < size(); ++i) {
            arrangeKidsByKeycoleAltitude(i);
        }
    }
    
    void fillKeyCols() {
        for (size_t i = 0; i < peaks.size(); ++i) {
            const auto &peak = peaks[i];
            if (peak.islandParent != kNoId && peak.islandParent != peak.id) {
                PeakId keyColId = extraPeakId(peak.keyColGeo());
                auto it = keyCols.insert({keyColId, {nodes[i].parent, {}}});
                assert(it.first->second.parent == nodes[i].parent);
                it.first->second.kids.push_back(i);
            }
        }
    }
    
    bool sameAltitudeParent(size_t i) {
        if (nodes[i].parent == kNoIndexSizeT || nodes[i].parent == i) {
            return false;
        }
        if (peaks[i].altitude == peaks[nodes[i].parent].altitude) {
            return true;
        }
        return false;
    }
    
    size_t size() const {
        assert(peaks.size() == nodes.size());
        return peaks.size();
    }
    
    void altitudeChangeConflicts(const std::map<PeakId, Altitude> &heightChange) {
        for (const auto &idHeight : heightChange) {
            auto it = idToIndex.find(idHeight.first);
            if (it == idToIndex.end()) {
                std::cout << idHeight.first << " not found\n";
                continue;
            }
            size_t index = it->second;
            const auto &peak = peaks[index];
            Altitude newHeight = idHeight.second;
            if (newHeight > peak.altitude) {
                continue;
            }
            size_t highestKeycolKid = kNoIndexSizeT;
            Altitude highestKeycol = kNoAltitude;
            for (size_t child : nodes[index].islandKids) {
                Altitude keycolAltitude = peaks[child].altitude - peaks[child].prominence;
                if (highestKeycolKid == kNoIndexSizeT || keycolAltitude > highestKeycol) {
                    highestKeycolKid = child;
                    highestKeycol = keycolAltitude;
                }
            }
            if (highestKeycolKid != kNoIndexSizeT) {
                if (highestKeycol > newHeight) {
                    std::cout << peak.toDescriptiveString() << "\n cant be lifted to " << newHeight << " because of (keycol: " << highestKeycol << ")\n" << peaks[highestKeycolKid].toDescriptiveString() << std::endl << std::endl;
                }
            }
            
        }
    }
    
    void switchParenthood(size_t parentIndex, size_t kidIndex) {
        arrangeKidsByKeycoleAltitude(parentIndex);
        arrangeKidsByKeycoleAltitude(kidIndex);
        auto &parent = peaks[parentIndex];
        auto &parentNode = nodes[parentIndex];
        auto kidIndexInParent = parentNode.kidIndexInParent(kidIndex);
        
        auto &kid = peaks[kidIndex];
        auto &kidNode = nodes[kidIndex];
        assert(parent.altitude <= kid.altitude);
        assert(kidIndexInParent != kNoIndexSizeT);
        if (parent.islandParent == parent.id) {
            kid.islandParent = kid.id;
            kidNode.parent = kidIndex;
        } else {
            kid.islandParent = parent.islandParent;
            kidNode.parent = parentNode.parent;
            if (!parent.isOrphan()) {
                auto parentIndexInParent = nodes[parentNode.parent].kidIndexInParent(parentIndex);
                assert(parentIndexInParent != kNoIndexSizeT);
                nodes[parentNode.parent].islandKids[parentIndexInParent] = kidIndex;
                arrangeKidsByKeycoleAltitude(parentNode.parent);
            }
        }
        
        parent.islandParent = kid.id;
        parentNode.parent = kidIndex;
        std::swap(kid.keyColLatitude, parent.keyColLatitude);
        std::swap(kid.keyColLongitude, parent.keyColLongitude);
        Altitude childKeycolAltitude = kid.altitude - kid.prominence;
        Altitude parentKeycolAltitude = parent.altitude - parent.prominence;
        
        assert(parentKeycolAltitude >= 0);
        kid.prominence = kid.altitude - parentKeycolAltitude;
        parent.prominence = parent.altitude - childKeycolAltitude;

        kidNode.islandKids.push_back(parentIndex);
        checkKidsOrder(kidIndex);
        kidNode.islandKids.insert(kidNode.islandKids.end(),
                                    parentNode.islandKids.begin() + kidIndexInParent + 1,
                                    parentNode.islandKids.end());
        parentNode.islandKids.erase(parentNode.islandKids.begin() + kidIndexInParent,
                                    parentNode.islandKids.end());
        checkKidsOrder(kidIndex);
        reIdKids(kidIndex);
        if (nodes[parentIndex].kidIndexInParent(kidIndex) != kNoIndexSizeT) {
            std::cout << parentIndex << endl;
            std::cout << kidIndex << endl;
            std::cout << parentNode.islandKids.size() << endl;
            for (const auto &kid : parentNode.islandKids) {
                std::cout << "kid: " << kid << endl;
            }
            assert(false);
        }
        for (size_t i = 1; i < kidNode.islandKids.size(); ++i) {
            if (peaks[kidNode.islandKids[i]].keyColAltitude() > peaks[kidNode.islandKids[i - 1]].keyColAltitude()) {
                assert(false);
            }
        }
        checkConsistency(parentIndex);
    }
    
    void adjustParenting(size_t index) {
        const auto &node = nodes.at(index);
        size_t parentIndex = node.parent;
        if (peaks.at(parentIndex).definitiveAltitudeLess(peaks.at(index))) {
            switchParenthood(parentIndex, index);
        }
    }
    
    void changeAltitude(PeakId id, Altitude newAltitude) {
        auto index = idToIndex[id];
        changeAltitudeByIndex(index, newAltitude);
    }
    
    void changeAltitudeByIndex(size_t index, Altitude newAltitude) {
        arrangeKidsByKeycoleAltitude(index);
        checkConsistency(index);
        Altitude oldAltitude = peaks[index].altitude;
        if (newAltitude == oldAltitude) {
            return;
        } else if (newAltitude > oldAltitude) {
            raiseAltitude(index, newAltitude);
        }
        else {
            flog() << "lowering " << peaks[index].name << " " << peaks[index].geo() << " alt from " << oldAltitude << " to " << newAltitude << std::endl;
            lowerAltitude(index, newAltitude);
        }
        arrangeKidsByKeycoleAltitude(index);
        checkConsistency(index);
    }
    
    void raiseAltitude(size_t index, Altitude newAltitude) {
        assert(peaks[index].altitude < newAltitude);
        Altitude keyColAltitude = peaks[index].keyColAltitude();
        peaks[index].altitude = newAltitude;
        peaks[index].prominence = newAltitude - keyColAltitude;
        for (;;) {
            auto &peak = peaks[index];
            if (peak.islandParent == peak.id || peak.islandParent == kNoId) {
                break;
            }
            auto &node = nodes[index];
            assert(node.parent != kNoIndexSizeT);
            auto &parent = peaks[node.parent];
            arrangeKidsByKeycoleAltitude(node.parent);
            auto &parentNode = nodes[node.parent];
            size_t indexInParent = parentNode.kidIndexInParent(index);
            assert(indexInParent != kNoIndexSizeT);
            if (peak.definitiveAltitudeLess(parent)) {
                break;
            }
            switchParenthood(node.parent, index);
        }
        checkConsistency(index);
    }
    
    bool checkProminence(size_t index) const {
        return peaks.at(index).altitude >= peaks.at(index).prominence;
    }
    
    void lowerAltitude(size_t index, Altitude newAltitude) {
        auto &node = nodes[index];
        auto &peak = peaks[index];
        assert(peak.altitude >= newAltitude);
        if (newAltitude < peak.keyColAltitude()) {
            std::cout << "cant lower: " << peak.toDescriptiveString() << " new_alt: " << newAltitude << " keycol_alt: " << peak.keyColAltitude() << endl;
            return;
        }
        if (node.islandKids.empty()) {
            Altitude keycolAlt = peak.keyColAltitude();
            peak.altitude = newAltitude;
            peak.prominence = newAltitude - keycolAlt;
            return;
        }
        arrangeKidsByKeycoleAltitude(index);
        while (newAltitude < peaks[node.islandKids.front()].keyColAltitude()
            && peaks[node.islandKids.front()].isExtra()) {
            rmTree(node.islandKids.front());
        }
        
        if (newAltitude < peaks[node.islandKids.front()].keyColAltitude()) {
            std::cout << "cant lift " << peak.id << " " << peak.name << " to " << newAltitude << std::endl;
            std::cout << peaks[node.islandKids.front()].keyColAltitude() << " min alt\n";
            std::cout << "problematic child: " << peaks[node.islandKids.front()].toDescriptiveString() << endl;
            return;
        }
        Altitude keyColAltitude = peak.keyColAltitude();
        peak.altitude = newAltitude;
        peak.prominence = newAltitude - keyColAltitude;
        
        size_t currentParent = index;
        std::vector<size_t> toCheckConsistency;
        toCheckConsistency.push_back(index);
        for (size_t kidsLeft = node.islandKids.size(); kidsLeft > 0; --kidsLeft) {
            auto &parentNode = nodes[currentParent];
            auto islandKidsCount = parentNode.islandKids.size();
            assert(islandKidsCount >= kidsLeft);
            auto &parent = peaks[currentParent];
            auto kidIndexInParent = islandKidsCount - kidsLeft;
            auto kidIndex = parentNode.islandKids[kidIndexInParent];
            assert(checkProminence(kidIndex));
            assert(checkProminence(currentParent));
            toCheckConsistency.push_back(kidIndex);
            auto &kid = peaks[kidIndex];

            assert(checkProminence(kidIndex));
            assert(checkProminence(currentParent));
            
            if (parent.altitude < kid.altitude || (parent.altitude == kid.altitude && parent.isExtra() && !kid.isExtra())) {
                switchParenthood(currentParent, kidIndex);
                checkConsistency(currentParent);
                currentParent = kidIndex;
            } else {
                checkConsistency(kidIndex);
            }
            assert(checkProminence(kidIndex));
            assert(checkProminence(currentParent));
        }
        checkConsistency(currentParent);
        for (size_t indexToCheck : toCheckConsistency) {
            checkConsistency(indexToCheck);
        }
    }
    
    void unnamePeak(size_t index) {
        auto &node = nodes[index];
        auto &peak = peaks[index];
        assert(!peak.isExtra());
        if (node.islandKids.empty()) {
            if (peak.id == peak.islandParent) {
                peak.islandParent = kNoId;
            }
            peak.id = kNoId;
            if (peak.islandParent == kNoId) {
                return;
            }
            auto parentNode = nodes[node.parent];
            bool deleted = false;
            for (auto it = parentNode.islandKids.begin();
                 it != parentNode.islandKids.end();
                 ++it) {
                if (*it == index) {
                    parentNode.islandKids.erase(it);
                    deleted = true;
                    break;
                }
            }
            assert(deleted);
            return;
        }
        PeakId newId = extraPeakId({peak.latitude, peak.longitude});
        if (idToIndex.count(newId)) {
            std::cout << "already have " << peaks.at(idToIndex.at(newId)).toDescriptiveString() << endl;
            std::cout << "current peaks: " << peak.toDescriptiveString() << std::endl;
//            peaks.id = kNoId;
            assert(false);
        }
        idToIndex.erase(peak.id);
        idToIndex[newId] = index;
        if (peak.islandParent == peak.id) {
            peak.islandParent = newId;
        }
        peak.id = newId;
        peak.name = "";
        reIdKids(index);
    }

    void arrangeKidsByKeycoleAltitude(size_t parent) {
        auto &node = nodes[parent];
        if (node.islandKids.empty()) {
            return;
        }
        std::sort(node.islandKids.begin(), node.islandKids.end(),
                  [&](size_t f, size_t s) {
            if (f == s) {
                return false;
            }
            return !peaks[f].definitiveKeyColAltitudeLess(peaks[s]);
        });
        checkKidsOrder(parent);
    }
    
    void reIdKids(size_t parent) {
        for (size_t kid : nodes.at(parent).islandKids) {
            peaks.at(kid).islandParent = peaks.at(parent).id;
            nodes.at(kid).parent = parent;
        }
    }
    
    void switchNames(size_t f, size_t s) {
        auto &fPeak = peaks[f];
        auto &sPeak = peaks[s];
        auto fId = fPeak.id;
        auto fName = sPeak.name;
        fPeak.id = sPeak.id;
        fPeak.name = sPeak.name;
        sPeak.id = fId;
        sPeak.name = fName;
        reIdKids(f);
        reIdKids(s);
    }
    
    void checkConsistency(size_t index) {
        if (peaks[index].id == kNoId) {
            return;
        }
        if (peaks[index].prominence > peaks[index].altitude) {
            throw std::logic_error("prominence more that altitude for " +
                peaks[index].toDescriptiveString());
        }
        if (nodes[index].parent != kNoIndexSizeT && nodes[index].parent != index) {
            if (peaks[index].altitude > peaks[nodes[index].parent].altitude) {
                throw std::logic_error(
                    peaks[index].toDescriptiveString() + "\n is higher than parent:\n" +
                    peaks[nodes[index].parent].toDescriptiveString());
            }
            if (nodes[nodes[index].parent].kidIndexInParent(index) == kNoIndexSizeT) {
                std::cout << "no kid index in parent for\n" << peaks[index].toDescriptiveString() << std::endl;
                throw std::logic_error("no kid index in parent for\n" +
                    peaks[index].toDescriptiveString());
            }
            assert(peaks[nodes[index].parent].id == peaks.at(index).islandParent);
        }
        if ((nodes[index].isOrphan()) != (peaks[index].isOrphan())) {
            std::cout << "node is orphan: " << nodes[index].isOrphan() << " peaks is orphan: " << peaks[index].isOrphan() << endl;
            throw std::logic_error(peaks.at(index).toDescriptiveString() + "\nnode is orphan: " + std::to_string(nodes[index].isOrphan())
                                   + " peaks is orphan: " + std::to_string(peaks[index].isOrphan()));
        }
    }
    
    void checkKidsOrder(size_t index) {
        assert(index < nodes.size());
        if (nodes[index].islandKids.empty()) {
            return;
        }
        for (auto it = nodes[index].islandKids.begin() + 1;
             it != nodes[index].islandKids.end(); ++it) {
            assert(*it != *(it - 1));
            if (peaks[*it].keyColAltitude() > peaks[*(it - 1)].keyColAltitude()) {
                std::cout << peaks[*it].toDescriptiveString() << std::endl;
                std::cout << peaks[*(it - 1)].toDescriptiveString() << std::endl;
                assert(false);
            }
        }
    }
    
    void checkConsistency() {
        assert(peaks.size() == nodes.size());
        for (size_t i = 0; i < peaks.size(); ++i) {
            checkConsistency(i);
        }
    }
    
    std::deque<size_t> sameAltitudeKids(size_t i) {
        std::deque<size_t> res;
        for (size_t k : nodes[i].islandKids) {
            assert(k != i);
            if (peaks[k].altitude == peaks[i].altitude) {
                res.push_back(k);
            }
        }
        return res;
    }
    
    void resolveSameAltitudeParenting(bool lowerExtras = true) {
        for (size_t index = 0; index < peaks.size(); ++index) {
            auto &node = nodes.at(index);
            if (node.parent == kNoIndexSizeT || node.parent == index) {
                continue;
            }
            
            auto &peak = peaks.at(index);
            auto &parentPeak = peaks.at(node.parent);
            if (parentPeak.definitiveAltitudeLess(peak, lowerExtras)) {
                switchParenthood(node.parent, index);
                assert(index > 0);
                index--;
            }
        }
        for (size_t index = 0; index < peaks.size(); ++index) {
            auto &node = nodes.at(index);
            if (node.parent == kNoIndexSizeT || node.parent == index) {
                continue;
            }
            auto &peak = peaks.at(index);
            auto &parentPeak = peaks.at(node.parent);
            assert(parentPeak.id == peak.islandParent);
            if (!peak.definitiveAltitudeLess(parentPeak, lowerExtras)) {
                std::cout << peak.toDescriptiveString() << endl;
                std::cout << "is higher than\n";
                std::cout << parentPeak.toDescriptiveString() << endl << endl;
            }
        }
    }
    
    void switchParentsIfNecessary(size_t index, size_t recursionDepth = 0) {
        if (recursionDepth > 10) {
            std::cout << "dbgrbd " << peaks[index].toDescriptiveString() << std::endl;
            assert(false);
        }
        if (nodes[index].parent == kNoIndexSizeT || nodes[index].parent == index) {
            return;
        }
        if (peaks[index].definitiveAltitudeLess(peaks[nodes[index].parent])) {
            return;
        }
        switchParenthood(nodes[index].parent, index);
        switchParentsIfNecessary(index, recursionDepth + 1);
    }
    
    void removeKid(size_t parent, size_t kid) {
        assert(parent != kid);
        auto &parentNode = nodes[parent];
        for (auto it = parentNode.islandKids.begin(); it != parentNode.islandKids.end(); ++it) {
            if (*it == kid) {
                nodes[kid].parent = kNoIndexSizeT;
                parentNode.islandKids.erase(it);
                return;
            }
        }
        std::cout << peaks[parent].toDescriptiveString() << endl;
        std::cout << peaks[kid].toDescriptiveString() << endl;
        assert(false);
    }
    
    void addKid(size_t parent, size_t kid) {
        auto &parentNode = nodes[parent];
        parentNode.islandKids.push_back(kid);
        nodes[kid].parent = parent;
        peaks[kid].islandParent = peaks[parent].id;
    }
    
    void resolveParentageInKeycols() {
        std::map<PeakId, std::vector<size_t>> keyColToPeaks;
        for (size_t i = 0; i < peaks.size(); ++i) {
            auto &peak = peaks[i];
            if (peak.islandParent == kNoId || peak.islandParent == peak.id) {
                continue;
            }
            PeakId colId = extraPeakId(peak.keyColGeo());
            auto it = keyColToPeaks.insert({colId, {}});
            it.first->second.push_back(i);
        }
        
        for (auto &keyColPeak : keyColToPeaks) {
            if (keyColPeak.second.size() < 2) {
                continue;
            }
            size_t parentIndex = nodes[keyColPeak.second.front()].parent;
            for (size_t i = 1; i < keyColPeak.second.size(); ++i) {
                size_t otherParentIndex = nodes[keyColPeak.second[i]].parent;
                if (otherParentIndex == parentIndex) {
                    continue;
                }
                if (peaks[parentIndex].definitiveAltitudeLess(peaks[otherParentIndex])) {
                    parentIndex = otherParentIndex;
                }
            }
            
            for (size_t i : keyColPeak.second) {
                size_t otherParentIndex = nodes[i].parent;
                if (otherParentIndex == parentIndex) {
                    continue;
                }
                removeKid(otherParentIndex, i);
                addKid(parentIndex, i);
            }
            arrangeKidsByKeycoleAltitude(parentIndex);
            checkConsistency(parentIndex);
        }
        
    }
    
    size_t nextNamedInTree(size_t index) const {
        std::vector<size_t> currentLevel;
        currentLevel.push_back(index);
        while (!currentLevel.empty()) {
            std::vector<size_t> nextLevel;
            for (size_t i : currentLevel) {
                for (size_t j : nodes[i].islandKids) {
                    if (!peaks[j].isExtra()) {
                        return j;
                    }
                    nextLevel.push_back(j);
                }
            }
            std::swap(currentLevel, nextLevel);
        }
        return kNoIndexSizeT;
    }
    
    Altitude highestCol(size_t index) {
        const auto &peak = peaks[index];
        Altitude res = peak.keyColAltitude();
        for (size_t kidIndex : nodes[index].islandKids) {
            const auto &kid = peaks[kidIndex];
            if (kid.keyColAltitude() > res) {
                res = kid.keyColAltitude();
            }
        }
        return res;
    }
    
    void rmTree(size_t index, bool rmNamed = false) {
        std::vector<size_t> currentLevel;
        currentLevel.push_back(index);
        int sanityCheck = 0;
        while (!currentLevel.empty()) {
            if (++sanityCheck > 1000) {
                throw std::logic_error("stuck in rm tree");
            }
            std::vector<size_t> nextLevel;
            for (size_t i : currentLevel) {
                if (!(rmNamed || peaks[i].isExtra())) {
                    std::cout << peaks[i].toDescriptiveString() << endl;
                    std::cout << peaks[nodes[i].parent].toDescriptiveString() << endl;
                    std::cout << peaks[index].toDescriptiveString() << endl;
                    throw std::logic_error("trying to rm named");
                }
                if (peaks[i].id == kNoId) {
                    throw std::logic_error("trying to rm already deleted");
                }
                peaks[i].id = kNoId;
                for (size_t j : nodes[i].islandKids) {
                    nextLevel.push_back(j);
                }
            }
            std::swap(currentLevel, nextLevel);
        }
        size_t parent = nodes.at(index).parent;
        if (parent != kNoIndexSizeT && parent != index) {
            removeKid(parent, index);
        }
    }
        
    void lowerUnnamed(size_t index) {
        const auto &peak = peaks[index];
        if (peak.isExtra()) {
            return;
        }
        static Altitude wl = 9000;
        if (peak.altitude < wl) {
            logTime("level");
            flog() << "lowering " << peak.id << " alt: " << peak.altitude << " geo: " << peak.geo() << endl;
            wl = peak.altitude;
        }
        auto &node = nodes[index];
        int sanityCheck = 0;
        bool fixedSomething = false;
        for (;;) {
            if (++sanityCheck > 10000) {
                std::cout << peak.toDescriptiveString() << "was raised 10k times\n";
                sanityCheck = 0;
            }
            if (node.parent == index || node.parent == kNoIndexSizeT) {
                break;
            }
            auto &parentPeak = peaks.at(node.parent);
            if (!parentPeak.isExtra()) {
                break;
            }
            size_t parentIndex = node.parent;
            lowerAltitude(parentIndex, peak.altitude);
            fixedSomething = true;
            if (parentIndex == node.parent) {
                throw std::logic_error("parentIndex == node.parent in lowerUnnamed");
            }
        }
        if (fixedSomething) {
            flog() << "fixed " << peak.toDescriptiveString() << endl;
        }
    }
    
    void lowerUnnamed(Altitude prominenceCutoff) {
        std::vector<size_t> indexes(peaks.size());
        for (size_t i = 0; i < peaks.size(); ++i) {
            indexes[i] = i;
        }
        
        std::sort(indexes.begin(), indexes.end(), [this](size_t f, size_t s) {
            return peaks[s].definitiveAltitudeLess(peaks[f]);
        });
        
        for (size_t i = 1; i < peaks.size(); ++i) {
            assert(peaks[indexes[i - 1]].altitude >= peaks[indexes[i]].altitude);
        }
        
        for (size_t i = 0; i < peaks.size(); ++i) {
            lowerUnnamed(indexes[i]);
        }
        for (size_t i = 0; i < peaks.size(); ++i) {
            if (peaks[i].isExtra()) {
                if (peaks[i].prominence < prominenceCutoff && peaks[i].id != kNoId) {
                    flog() << "rm " << peaks[i].id << std::endl;
                    rmTree(i);
                    flog() << "rm " << peaks[i].id << " finished" << std::endl;
                }
                continue;
            }
            if (nodes[i].parent != kNoIndexSizeT && peaks[nodes[i].parent].isExtra()) {
                throw std::logic_error(peaks[i].toDescriptiveString() + " has extraPeak parent\n");
            }
        }
    }

    void correctIsolationByNHN() {
        for (auto &peak : peaks) {
            if (peak.nhn == kNoId) {
                continue;
            }
            const auto &nhn = peaks.at(idToIndex.at(peak.nhn));
            if (peak.ilpLatitude == kNoDegrees ||
                DistanceTools::distanceOnEarth(peak.geo(), nhn.geo()) < DistanceTools::distanceOnEarth(peak.geo(), peak.ilp()))
            {
                peak.ilpLatitude = nhn.latitude;
                peak.ilpLongitude = nhn.longitude;
            }
        }
    }
    
};
