#pragma once

#include "index.h"
#include "dem_grid.h"
#include "ridge.h"
#include "bitfield.h"

#include <unordered_map>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include <deque>
#include <functional>
#include <memory>
#include <fstream>

//#define SIMPLE_RIDGE_CALCULATOR_DEßBUG
#ifdef SIMPLE_RIDGE_CALCULATOR_DEBUG
#define srica_dbg(NX) NX
#else
#define srica_dbg(NX)
#endif

using AltitudeGetter = std::function<Altitude(const PlanarIndex&)>;

template <typename View>
struct RidgeCalculator {
    PlanarIndexDir kPeakGradDirection = kPlanarIndexUp;
    static constexpr Altitude kProminenceCutoff = 10;

    struct RidgeInfo {
        Ridge ridgeToFirst;
        uint32_t parentNodeIndex = k32NoIndex;
        uint32_t firstNodeIndex = k32NoIndex;
        PlanarIndexDir dirFromCol = kPlanarIndexNoDir;
        uint8_t edgeMask = 0;

//        RidgeInfo(RidgeInfo &&other) : ridgeToFirst(std::move(other.ridgeToFirst)), parentNodeIndex(other.parentNodeIndex),
//                                       firstNodeIndex(other.firstNodeIndex), dirFromCol(other.dirFromCol), edgeMask(other.edgeMask){}

    };

    struct PeakNode {
        const SpatialIndex spatial;
        uint32_t islandParentIndex = k32NoIndex;
        uint32_t keyColIndex = k32NoIndex;
        uint32_t nextInLineCol = k32NoIndex;
        uint32_t peakTreeIndex = k32NoIndex;

        std::vector<uint32_t> indexRidge;
        uint8_t edgeMask = 0;
        uint32_t ridgesCount = 0;
    };

    struct ColNode {
        const SpatialIndex spatial;
        std::vector<uint32_t> adjacentPeaksIndexes;
        std::vector<Ridge> ridges;
        std::vector<PlanarIndexDir> peakDirs;
        uint32_t nextInLineIndexInCol = k32NoIndex;
        uint8_t edgeMask = 0;

        void addRidge(RidgeInfo &ridgeInfo) {
            uint32_t peakIndex = ridgeInfo.firstNodeIndex;
            adjacentPeaksIndexes.push_back(peakIndex);
            ridgeInfo.ridgeToFirst.checkConsistency();
            ridges.emplace_back(std::move(ridgeInfo.ridgeToFirst));
            ridges.back().checkConsistency();
            peakDirs.push_back(ridgeInfo.dirFromCol);
        }

        [[nodiscard]] bool empty() const {
            return adjacentPeaksIndexes.empty();
        }

        [[nodiscard]] uint32_t nextInLineNodeIndex() const {
            assert(nextInLineIndexInCol != k32NoIndex);
            return adjacentPeaksIndexes.at(nextInLineIndexInCol);
        }

        void setNextInLineNode(uint32_t nodeIndex) {
            nextInLineIndexInCol = peakIndexInCol(nodeIndex);
        }

        [[nodiscard]] const Ridge& ridgeToPeak(uint32_t peakIndex) const {
            return ridges.at(peakIndexInCol(peakIndex));
        }

        [[nodiscard]] uint8_t peakIndexInCol(uint32_t peakIndex) const {
            for (uint8_t i = 0; i < adjacentPeaksIndexes.size(); ++i) {
                if (adjacentPeaksIndexes.at(i) == peakIndex) {
                    return i;
                }
            }
            assert(false);
        }
    };

    struct PotentialCol {
        SpatialIndex spatial;
        uint8_t directionsUpMask;
        uint8_t edgeMask;
    };

    struct GradientWalkResult {
        PlanarIndex endPlanar;
        uint8_t edgeMask = 0;
    };

    uint64_t totalRidgesLength = 0;
    uint64_t totalPeaksChecked = 0;

    std::vector<PeakNode> peakNodes;
    std::vector<ColNode> colNodes;
    std::vector<PotentialCol> potentialCols;
    CompactGrid<View, 1> grad;
    BoolDequeField peakField;
    std::unordered_map<PlanarIndex, uint32_t, PlanarIndexHash> planarToNodeIndex;
    std::unique_ptr<AltitudeGetter> altitudeGetter = nullptr;
    std::deque<Peak> *peaksp;
    Converter *converterp;
    Rectangle innerArea;
    const std::string ridgesDir;

    template <typename DEMView>
    RidgeCalculator(const DEMGrid<DEMView> &region, std::string ridgesDir, const Index &margin = 1)
        : grad("", region.cellsFrom - region.positivityShift, region.cellsTo - region.positivityShift, region.cellsResolution, region.cellsResolution, false, true),
        peakField(region.indexFrom(), region.indexTo()),
        ridgesDir(std::move(ridgesDir))
    {
        innerArea = Rectangle{region.indexFrom() + region.cellsResolution * margin,
            region.indexTo() - region.cellsResolution * margin};
        std::cout << "inner origin: " << innerArea.origin << " inner end: " << innerArea.end << endl;
        altitudeGetter = std::make_unique<AltitudeGetter>([&region](const PlanarIndex &planar) {
            return region.get(planar);
        });

        std::cout << "dbg planar: " << debugPlanar << endl;

        uint32_t peaksCounter = 0;

        region.doForEachPlanar([&region, &peaksCounter, this](const PlanarIndex &planar) {
            SpatialIndex spatial{planar, region.get(planar)};
            if (spatial.altitude <= region.kFallBackValue()) {
                return;
            }

            if (spatial.planar == debugPlanar) {
                std::cout << "debug spatial: " << spatial << endl;
            }

            auto isPotentiallyCol = region.template isPotentialCol<1>(spatial, true);
            auto isPeak = region.runningVicinityInfo.highestNeighborDir == kPlanarIndexNoDir;
            assert(isPeak == (region.runningVicinityInfo.dirsUpMask == 0));

            assert(!(isPeak && isPotentiallyCol));

            uint8_t gradDirection = isPeak ? kPeakGradDirection : region.runningVicinityInfo.highestNeighborDir;

            if (isPeak) {
                peakField.set(planar, true);
                ++peaksCounter;
            }

            grad.set(planar, gradDirection);
            assert(grad.get(planar) == gradDirection);

            if (planar == debugPlanar) {
                std::cout << "dbgDir: " << gradDirection << endl;
                std::cout << "grad at dbg planar: " << grad.get(debugPlanar) << endl;
            }

            if (region.runningVicinityInfo.edgeMask) {
                if (!isPeak) {
                    if ((offsetMask(PlanarIndexDir(gradDirection)) & region.runningVicinityInfo.edgeMask) != 0) {
                        std::cout << spatial << " " << region.planarToGeodesic(planar) << endl;
                        assert(false && "facing out");
                    }
                    potentialCols.push_back({spatial, region.runningVicinityInfo.dirsUpMask, region.runningVicinityInfo.edgeMask});
                }
                return;
            }

            if (isPotentiallyCol) {
                if (region.template isPotentialCol<3>(spatial)) {
                    potentialCols.push_back({spatial, region.runningVicinityInfo.dirsUpMask});
                }
            }
        });
        std::cout << "total cols: " << potentialCols.size() << endl;
        std::cout << "total peaks: " << peaksCounter << endl;
    }

    void addPeaksTree(std::deque<Peak> &peaks, Converter &converter) {
        peaksp = &peaks;
        converterp = &converter;
        uint32_t addedPeaks = 0;
        for (uint32_t i = 0; i < peaks.size(); ++i) {
            const auto &peak = peaks[i];

            if (peak.prominence == 0) {
                continue;
            }

            PlanarIndex peakPlanar = converter.geodesicToPlanar(peak.geo());
            uint32_t nodeIndex = addNode({peakPlanar, peak.altitude});
            if (nodeIndex != k32NoIndex) {
                peakNodes.at(nodeIndex).peakTreeIndex = i;
                ++addedPeaks;
            }

        }
        std::cout << "added " << addedPeaks << " peaks\n";
    }

    void checkField() {

    }

    void calculate() {
        std::sort(potentialCols.rbegin(), potentialCols.rend(), [](const PotentialCol &f, const PotentialCol &s) {
            return f.spatial < s.spatial;
        });

        for (size_t i = 1; i < potentialCols.size(); ++i) {
            assert(potentialCols[i].spatial < potentialCols[i - 1].spatial);
        }
        Altitude waterLevel = kHighNansAltitude;
        for (uint32_t colIndex = 0; colIndex < potentialCols.size(); ++colIndex) {
            auto &col = potentialCols[colIndex];
            if (col.spatial.altitude < waterLevel) {
                waterLevel = col.spatial.altitude;
                std::cout << "wl: " << waterLevel << endl;
            }
            calculateCol(col);
        }

        uint32_t totalPeaks = 0;
        uint32_t smallPeaks = 0;
        uint32_t totalCols = 0;
        uint32_t isolatedPeaks = 0;
        uint32_t isolatedSmallPeaks = 0;
        uint32_t triPeaks = 0;
        uint32_t smallTriPeaks = 0;
        uint32_t fatPeaks = 0;
        uint32_t smallFatPeaks = 0;

        for (const auto &peak : peakNodes) {
            if (peak.spatial.planar.isInsideRectangle(innerArea)) {
                totalPeaks++;
                if (getProminence(peak) < 100) {
                    smallPeaks++;
                }
                if (peak.ridgesCount < 3) {
                    if (getProminence(peak) < 100) {
                        isolatedSmallPeaks++;
                    }
                    isolatedPeaks++;
                } else if (peak.ridgesCount == 3) {
                    if (getProminence(peak) < 100) {
                        ++smallTriPeaks;
                    }
                    triPeaks++;
                } else {
                    if (getProminence(peak) < 100) {
                        ++smallFatPeaks;
                    }
                    fatPeaks++;
                }
            }
        }
        for (const auto &col : colNodes) {
            if (col.spatial.planar.isInsideRectangle(innerArea)) {
                totalCols++;
            }
        }

        std::cout << OUTSTR(totalRidgesLength) << std::endl;
        std::cout << OUTSTR(totalPeaksChecked) << endl;
        std::cout << "totalPeaks:\t" << totalPeaks << "\nsmallPeaks:\t" << smallPeaks << "\ntotalCols:\t" << totalCols << "\nisolatedPeaks:\t" << isolatedPeaks << "\ntriPeaks:\t" << triPeaks << "\nsmallTriPeaks:\t" << smallTriPeaks << "\nfatPeaks:\t" << fatPeaks << "\nsmallFatPeaks:\t" << smallFatPeaks << endl;
        std::cout << "real peaks left: " <<  (totalPeaks - isolatedSmallPeaks) << endl;
    }

    void calculateCol(PotentialCol &col) {
        srica_dbg(std::cout << "calculating col: " << col.spatial << endl;)
        const auto &colPlanar = col.spatial.planar;
        char colEdgeMask = col.edgeMask;

        std::map<SpatialIndex, RidgeInfo> ridgeInfos;

        for (uint8_t offsetIndex = 0; offsetIndex < kMaxNeighborsCount; ++offsetIndex) {
            auto offsetDir = PlanarIndexDir(offsetIndex);
            if (!checkMask(col.directionsUpMask, offsetIndex)) {
                continue;
            }
            calculateRidgeFromColInDir(ridgeInfos, colPlanar, offsetDir);
        }

        for (const auto &spatialRidgeInfo : ridgeInfos) {
            colEdgeMask |= spatialRidgeInfo.second.edgeMask;
            spatialRidgeInfo.second.ridgeToFirst.checkConsistency();
            assert(spatialRidgeInfo.second.parentNodeIndex != k32NoIndex || spatialRidgeInfo.second.edgeMask == 0);
        }

        assert(!ridgeInfos.empty());
        auto ridgeIt = ridgeInfos.begin();

        while (ridgeIt != std::prev(ridgeInfos.end())) {
            auto &ridgeInfo = ridgeIt->second;
            const auto &ridgeParentSpatial = ridgeIt->first;
            Altitude prominence = ridgeParentSpatial.altitude - col.spatial.altitude;
            if (prominence < kProminenceCutoff && ridgeInfo.parentNodeIndex == k32NoIndex) {
                reverseRidge(col.spatial.planar, ridgeParentSpatial.planar, ridgeInfo.dirFromCol);
                ridgeIt = ridgeInfos.erase(ridgeIt);
            } else {
                ++ridgeIt;
            }
        }

        for (const auto &spatialRidgeInfo : ridgeInfos) {
            spatialRidgeInfo.second.ridgeToFirst.checkConsistency();
        }

        auto parentIt = std::prev(ridgeInfos.end());
        const SpatialIndex &parentSpatial = parentIt->first;
        RidgeInfo &parentInfo = parentIt->second;

        for (const auto &spatialRidgeInfo : ridgeInfos) {
            spatialRidgeInfo.second.ridgeToFirst.checkConsistency();
        }

        grad.set(colPlanar, parentInfo.dirFromCol);

        if (ridgeInfos.size() == 1) {
            if (colEdgeMask != 0) {
                if (parentInfo.parentNodeIndex == k32NoIndex) {
                    parentInfo.parentNodeIndex = addNode(parentSpatial);
                }
                peakNodes[parentInfo.parentNodeIndex].edgeMask |= colEdgeMask;
            }
            srica_dbg(std::cout << col.spatial << " is not a keyCol\n";)
            return;
        }

        uint32_t colNodeIndex = addColNode(col.spatial);
        auto &colNode = colNodes.at(colNodeIndex);

        for (const auto &spatialRidgeInfo : ridgeInfos) {
            spatialRidgeInfo.second.ridgeToFirst.checkConsistency();
        }

        for (auto &spatialRidgeInfo : ridgeInfos) {
            auto &ridgeInfo = spatialRidgeInfo.second;
            if (ridgeInfo.parentNodeIndex == k32NoIndex) {
                ridgeInfo.parentNodeIndex = addNode(spatialRidgeInfo.first);
                ridgeInfo.firstNodeIndex = ridgeInfo.parentNodeIndex;
            }
            peakNodes.at(ridgeInfo.parentNodeIndex).edgeMask = colEdgeMask;
            peakNodes.at(ridgeInfo.parentNodeIndex).ridgesCount++;
            ridgeInfo.ridgeToFirst.checkConsistency();
            colNode.addRidge(ridgeInfo);
        }
        for (const auto &r : colNode.ridges) {
            r.checkConsistency();
        }

        uint32_t parentIndex = parentInfo.parentNodeIndex;
        assert(parentIndex != k32NoIndex);
        colNode.setNextInLineNode(parentInfo.firstNodeIndex);

        std::vector<uint32_t> parentIndexRidge = getIndexRidge(colNodeIndex, parentInfo.firstNodeIndex, parentInfo.parentNodeIndex, false);

//        srica_dbg(std::cout << "parent index ridge:\n"; logIndexRidge(parentIndexRidge, true);)
        for (const Ridge &ridge : colNode.ridges) {
            totalRidgesLength += ridge.length();
        }

        for (auto it = ridgeInfos.begin(); it != std::prev(ridgeInfos.end()); ++it) {
            auto subductedNodeIndex = it->second.parentNodeIndex;
            if (subductedNodeIndex == k32NoIndex) {
                assert(false);
            }
            auto &subductedNode = peakNodes.at(subductedNodeIndex);
            srica_dbg(std::cout << "subducting " << subductedNodeIndex << " " << subductedNode.spatial << " bride: " << it->second.firstNodeIndex << " parent: " << parentIndex << " " << peakNodes.at(parentIndex).spatial << endl;)
            auto &subductedIndexRidge = subductedNode.indexRidge;
            subductedIndexRidge = getIndexRidge(colNodeIndex, it->second.firstNodeIndex, it->second.parentNodeIndex, true);
//            srica_dbg(std::cout << "subductedIndexRidge to col\n"; logIndexRidge(subductedIndexRidge, true));
            subductedIndexRidge.push_back(colNodeIndex);
            assert(subductedIndexRidge.size() % 2 == 0);
            subductedIndexRidge.insert(subductedIndexRidge.end(), parentIndexRidge.begin(), parentIndexRidge.end());
            subductedNode.islandParentIndex = parentIndex;
            subductedNode.keyColIndex = colNodeIndex;
//            if (subductedNode.spatial.planar.isInsideRectangle(safeArea)) {
//                std::cout << subductedNode.spatial << " " << colNode.spatial << " " << getProminence(subductedNode) << endl;
//            }
//            srica_dbg(std::cout << "subductedNode.IndexRidge\n"; logIndexRidge(subductedIndexRidge, true));
            assert(subductedNode.nextInLineCol != k32NoIndex);
            if (subductedNode.peakTreeIndex != k32NoIndex && subductedNode.edgeMask == 0) {
                serializePeakNodeRidge(subductedNodeIndex);
                auto &peak = peaksp->at(subductedNode.peakTreeIndex);
                if (getProminence(subductedNode) != peak.prominence) {
                    std::cout << "inconsistent node\n";
                    std::cout << peak.toDescriptiveString() << endl;
                    std::cout << "node: " << subductedNode.spatial <<  " col: " << colNode.spatial << " parent: " << peakNodes.at(parentIndex).spatial << endl;
                } else {
                    ++totalPeaksChecked;
                }
            }
        }
    }

    Altitude getProminence(const PeakNode &node) const {
        if (node.keyColIndex == k32NoIndex) {
            return kNoAltitude;
        }
        const auto &colAlt = colNodes.at(node.keyColIndex).spatial.altitude;
        assert(colAlt > 0);
        return node.spatial.altitude - colAlt;
    }

    std::vector<uint32_t> getIndexRidge(uint32_t colIndex, uint32_t peakIndex,
                                        uint32_t parentPeakIndex, bool reverseWhileGetting) {
        assert(parentPeakIndex != k32NoIndex && colIndex != k32NoIndex && peakIndex != k32NoIndex);
        srica_dbg(std::cout << "getIndexRidge col " << colIndex << " peaks " << peakIndex << " parentIndex: " << parentPeakIndex << endl;)
        std::vector<uint32_t> indexRidge;
        indexRidge.push_back(peakIndex);

        while (peakIndex != parentPeakIndex) {
            assert(peakNodes[peakIndex].nextInLineCol != k32NoIndex);
            assert(peakNodes.at(peakIndex).spatial < peakNodes.at(parentPeakIndex).spatial);
            uint32_t nextColIndex = peakNodes.at(peakIndex).nextInLineCol;
            uint32_t nextNodeIndex = colNodes.at(nextColIndex).nextInLineNodeIndex();

            if (reverseWhileGetting) {
                peakNodes.at(peakIndex).nextInLineCol = colIndex;
                colNodes.at(nextColIndex).setNextInLineNode(peakIndex);
            }
            peakIndex = nextNodeIndex;
            colIndex = nextColIndex;
            indexRidge.push_back(colIndex);
            indexRidge.push_back(peakIndex);
        }

        if (reverseWhileGetting) {
            peakNodes.at(peakIndex).nextInLineCol = colIndex;
            std::reverse(indexRidge.begin(), indexRidge.end());
        }

        return indexRidge;
    }

    void logIndexRidge(const std::vector<uint32_t> &indexRidge, bool startFromPeak) const {
        bool isPeak = startFromPeak;
        for (const auto &i : indexRidge) {
            std::cout << i << (isPeak ? " peaks " : " col  ") << (isPeak ? peakNodes.at(i).spatial : colNodes.at(i).spatial)
                      << " next: " << (isPeak ? peakNodes.at(i).nextInLineCol : colNodes.at(i).nextInLineNodeIndex()) << endl;
            isPeak = !isPeak;
        }
        std::cout << endl;
    }

    void calculateRidgeFromColInDir(std::map<SpatialIndex, RidgeInfo> &ridgeInfos, const PlanarIndex &colPlanar, PlanarIndexDir offsetDir) {
        if (colPlanar == debugPlanar) {
            std::cout << "going in dir " << offsetDir << endl;
        }
        Ridge ridge(colPlanar);
        auto walkResult = gradientWalk(colPlanar, offsetDir, [&ridge, this](const PlanarIndex& currentPlanar, PlanarIndexDir dir, PlanarIndexDir previousDir) {
            ridge.addStep(currentPlanar);
            return !peakField.get(currentPlanar);

        });
        const auto &peakPlanar = walkResult.endPlanar;

        auto findIt = planarToNodeIndex.find(peakPlanar);
        const uint32_t firstNodeIndex = findIt == planarToNodeIndex.end() ? k32NoIndex : findIt->second;
        uint32_t nodeIndex = firstNodeIndex;
        if (nodeIndex != k32NoIndex) {
            uint32_t dbgCounter = 0;
            for (;;) {
                if (++dbgCounter > 100000) {
                    assert(false);
                }
                if (peakNodes[nodeIndex].islandParentIndex != k32NoIndex) {
//                    walkResult.edgeMask |= peakNodes[nodeIndex].edgeMask;
                    nodeIndex = peakNodes[nodeIndex].islandParentIndex;
                } else {
                    break;
                }
            }
        }

        SpatialIndex peakSpatial;

        if (nodeIndex == k32NoIndex) {
            assert(altitudeGetter != nullptr);
            peakSpatial = SpatialIndex{peakPlanar, (*altitudeGetter)(peakPlanar)};
        } else {
            peakSpatial = peakNodes[nodeIndex].spatial;
        }

        if (colPlanar == debugPlanar) {
            std::cout << "going in dir " << offsetDir << " up to " << peakSpatial << endl;
        }

        uint8_t edgeMask = nodeIndex == k32NoIndex ? 0 : peakNodes.at(nodeIndex).edgeMask;
        ridge.checkConsistency();
        std::pair<SpatialIndex, RidgeInfo> pairToInsert{peakSpatial, {std::move(ridge), nodeIndex, firstNodeIndex, offsetDir, edgeMask}};
        pairToInsert.second.ridgeToFirst.checkConsistency();
        auto ins = ridgeInfos.insert(pairToInsert);
        if (!ins.second) {
            if (pairToInsert.second.ridgeToFirst.isShorter(ins.first->second.ridgeToFirst)) {
                ins.first->second = pairToInsert.second;
            }
        }
    }

    GradientWalkResult gradientWalk(const PlanarIndex &start, PlanarIndexDir startingDir,
                                    const std::function<bool(const PlanarIndex&, PlanarIndexDir, PlanarIndexDir)> &stepCallback,
                                    bool isRidgeWalk = true) {
        PlanarIndex previousPlanar = start;
        PlanarIndexDir previousDir = startingDir;
        PlanarIndex currentPlanar = nextInDir(start, startingDir);
        PlanarIndexDir currentDir = PlanarIndexDir(grad.get(currentPlanar));
        GradientWalkResult result{start, grad.getEdgeMask(start)};

        auto startSpatial = SpatialIndex{start, (*altitudeGetter)(start)};

        uint32_t dbgCounter = 0;
        for (;;) {
            if (++dbgCounter > 100000) {
                std::cout << "looped\n";
                assert(false);
            }

//            char edgeMask = grad.getEdgeMask(currentPlanar);
//            if (edgeMask) {
//                result.edgeMask |= edgeMask;
//            }

//            SpatialIndex currentSpatial{currentPlanar, (*altitudeGetter)(currentPlanar)};
//            if (isRidgeWalk && !(startSpatial < currentSpatial)) {
//                std::cout << currentSpatial << " currentDir: " << currentDir << " lower than " << startSpatial << " at: " << start << endl;
//                std::cout << "previous: " << previousPlanar << " dir: " << previousDir << endl;
//                assert(false && "lower next alt");
//            }

//            if (currentPlanar == debugPlanar) {
//
//            }
            if (!stepCallback(currentPlanar, currentDir, previousDir)) {
                result.endPlanar = currentPlanar;
                break;
            };

            previousPlanar = currentPlanar;
            previousDir = currentDir;

            currentPlanar = nextInDir(previousPlanar, previousDir);

//            if (edgeMask && !grad.planarIsInside(currentPlanar)) {
//                result.endPlanar = previousPlanar;
//                if (previousDir != kPeakGradDirection) {
//                    std::cout << previousPlanar << " " << previousDir << " from " << start << endl;
//                    assert(false);
//                }
//                break;
//            }
            assert(grad.planarIsInside(currentPlanar));

            currentDir = PlanarIndexDir(grad.get(currentPlanar));
        }
        assert(result.endPlanar != start);

        if (isRidgeWalk) {
            SpatialIndex endSpatial{result.endPlanar, (*altitudeGetter)(result.endPlanar)};
            assert(startSpatial < endSpatial);
        }
        return result;
    }

    uint32_t addNode(const SpatialIndex &spatial, char edgeMask = 0) {
        uint32_t nodeIndex = static_cast<uint32_t>(peakNodes.size());
        auto ins = planarToNodeIndex.insert({spatial.planar, nodeIndex});
        if (!ins.second) {
            return k32NoIndex;
        }
        peakNodes.push_back(PeakNode{spatial});
        peakNodes.back().edgeMask = edgeMask;
        return nodeIndex;
    }

    uint32_t addColNode(const SpatialIndex &spatial) {
        uint32_t nodeIndex = static_cast<uint32_t>(colNodes.size());
        colNodes.emplace_back(ColNode{spatial});
        srica_dbg(std::cout << "adding col " << spatial << " at index " << nodeIndex << endl;)
        return nodeIndex;
    }

    void reverseRidge(const PlanarIndex &from, const PlanarIndex &to, PlanarIndexDir dir) {
//        std::cout << "reversing ridge from " << from << " to " << to << endl;
        auto walkResult = gradientWalk(from, dir, [this, &to](const PlanarIndex &currentPlanar, PlanarIndexDir currentDir, PlanarIndexDir previousDir) {
            grad.set(currentPlanar, previousDir ^ 1);
            return currentPlanar != to;
        });

        gradientWalk(to, PlanarIndexDir(grad.get(to)), [&from](const PlanarIndex &currentPlanar, PlanarIndexDir currentDir, PlanarIndexDir previousDir) {
            return currentPlanar != from;
        }, false);

        peakField.set(to, false);

        assert(walkResult.endPlanar == to);
    }

    int serializePeakNodeRidge(uint32_t nodeIndex) {
        const auto &node = peakNodes.at(nodeIndex);
        if (node.islandParentIndex == k32NoIndex) {
            return -1;
        }

        if (node.peakTreeIndex == k32NoIndex) {
            return -2;
        }

        if (node.edgeMask) {
            return -3;
        }

        const auto &peak = peaksp->at(node.peakTreeIndex);
        std::string filename = ridgesDir + std::to_string(peak.id) + ".txt";
        std::ofstream out(filename);
//        std::cout << "serializing " << filename << endl;
        const auto &indexRidge = peakNodes.at(nodeIndex).indexRidge;
        uint32_t lastIndexIndex = static_cast<uint32_t>(indexRidge.size() - 1);
        srica_dbg(std::cout << "serializing ridge from " << node.spatial << " to " << peakNodes.at(indexRidge.back()).spatial << " index ridge length " << indexRidge.size() << endl;);
//        outputPlanar(out, node.spatial.planar, node.peakTreeIndex);
        for (uint32_t i = 0; i < lastIndexIndex; i += 2) {
            const auto &peakNode = peakNodes.at(indexRidge.at(i));
            outputPlanar(out, peakNode.spatial.planar, peakNode.peakTreeIndex);
            if (peakNodes.at(indexRidge.at(i)).spatial.altitude > node.spatial.altitude) {
                srica_dbg(std::cout << "got line parent\n");
                return 0;
            }

            const auto &ridgeToCol = ridgeBetweenNodes(indexRidge.at(i + 1), indexRidge.at(i));
            srica_dbg(std::cout << "ridge to col size: " << ridgeToCol.length() << endl;);
            ridgeToCol.testWalkThrough();
            ridgeToCol.walkThroughEndsExclusive([this, &out](const PlanarIndex &planar) {
                outputPlanar(out, planar);
            }, true);
            const auto &ridgeFromCol = ridgeBetweenNodes(indexRidge.at(i + 1), indexRidge.at(i + 2));
            ridgeFromCol.checkConsistency();
            ridgeFromCol.testWalkThrough();
            const auto &colNode = colNodes.at(indexRidge.at(i + 1));
            outputPlanar(out, colNode.spatial.planar, k32NoIndex);
            ridgeFromCol.walkThroughEndsExclusive([this, &out](const PlanarIndex &planar) {
                outputPlanar(out, planar);
            }, false);
        }
        auto &parentNode = peakNodes.at(indexRidge.back());
        outputPlanar(out, parentNode.spatial.planar, parentNode.peakTreeIndex);
        return 0;
    }

    void outputPlanar(std::ofstream &out, const PlanarIndex &planar, uint32_t peakTreeIndex = k32NoIndex, bool markCol = false) {
        auto geo = converterp->planarToGeodesic(planar);
        out << geo.x << " " << geo.y << " " << (*altitudeGetter)(planar);
        if (peakTreeIndex != k32NoIndex) {
            auto &peak = peaksp->at(peakTreeIndex);
            out << " " << peak.id << " " << peak.altitude << " " << peak.prominence << " " << peak.name;
        } else if (markCol) {
            out << " " << "col";
        }
        out << endl;
        if (!out) {
            std::cout << "broken out in outputPlanar\n";
            assert(false);
        }
    }

    const Ridge& ridgeBetweenNodes(uint32_t colIndex, uint32_t peakIndex) const {
        return colNodes.at(colIndex).ridgeToPeak(peakIndex);
    }
//    struct ColNode {
//            const SpatialIndex spatial;
//            std::vector<uint32_t> adjacentPeaksIndexes;
//            std::vector<Ridge> ridges;
//            std::vector<PlanarIndexDir> peakDirs;
//            uint32_t nextInLineIndexInCol = k32NoIndex;
//            char edgeMask = 0;
//
//            void addRidge(RidgeInfo &ridgeInfo);
//            bool empty() const;
//            uint32_t nextInLineNodeIndex();
//            void setNextInLineNode(uint32_t nodeIndex);
//            const Ridge& ridgeToPeak(uint32_t peakIndex);
//            uint8_t peakIndexInCol(uint32_t peakIndex) const;
//    };
    void prune(Altitude prominenceCutoff) {

        for (int i = 0; i < 5; ++i) {
            logTime();

            for (Index colIndex = 0; colIndex < colNodes.size(); ++colIndex) {

            }

            logTime("pruning iteration");
        }
    }

};
