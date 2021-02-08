#pragma once

#include "index.h"
#include "mapped_file.h"
#include "field_view.h"
#include "dem_filename_formatter.h"
#include "bitfield.h"
#include "../Lib/converter.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <climits>
#include <sstream>

template <typename View = FieldView, typename Value = Altitude>
struct PlanarGrid {
    struct Vicinity {
        std::vector<SpatialIndex> neighbors;
        uint16_t highestNeighbor = 0;
        uint16_t lowestNeighbor = 0;
        uint16_t edgeMask = 0;
    };

    const PlanarIndex cellsResolution;
    PlanarIndex cellsFrom, cellsTo;
    PlanarIndex viewsSize;
    const PlanarIndex positivityShift;
    const PlanarIndex fileSize;
    std::vector<MappedFile> files;
    std::vector<View> views;
    std::vector<std::vector<Value>> buffers;
    unsigned filesCount = 0;
    PlanarIndex _indexFrom;
    PlanarIndex _indexTo;
    Converter converter;
    bool verbose;

    PlanarGrid(const std::string &directoryPath,
               const PlanarIndex &from = kLatLngFrom,
               const PlanarIndex &to = kLatLngTo,
               PlanarIndex fileSize = {1200, 1200},
               PlanarIndex cellResolution = {1200, 1200},
               bool verbose = false,
               bool storeHeightfield = false)
    : cellsResolution(cellResolution),
    positivityShift{-from.x, -from.y},
    fileSize(fileSize),
    converter(positivityShift, cellsResolution),
    verbose(verbose) {
        init(directoryPath, from, to, storeHeightfield);
    }

    PlanarGrid(const std::string &directoryPath,
                const PlanarIndex &from = kLatLngFrom,
                const PlanarIndex &to = kLatLngTo,
                PlanarIndex fileSize = {1200, 1200},
                Index cellResolution = 1200,
                bool verbose = false,
                bool storeHeightfield = false)
    : PlanarGrid(directoryPath, from, to, fileSize, {cellResolution, cellResolution}, verbose, storeHeightfield) {}

    ~PlanarGrid() {
        deinit();
    }

    void init(const std::string &directoryPath, const PlanarIndex &from, const PlanarIndex &to, bool storeHeightfield) {
        deinit();
        if (verbose) {
            std::cout << "creating dem grid dir: " << directoryPath << std::endl;
        }
        bool creatingNew = directoryPath.empty();
        if (creatingNew) assert(storeHeightfield);

        cellsFrom = from + positivityShift;
        cellsTo = to + positivityShift;

        if (cellsFrom != PlanarIndex{0, 0}) {
            assert(false);
        }

        filesCount = 0;
        uint16_t missingFiles = 0;

        viewsSize = PlanarIndex{cellsTo.latitude - cellsFrom.latitude + 1, cellsTo.longitude - cellsFrom.longitude + 1};
        size_t viewsCount = viewsSize.x * viewsSize.y;
        views.resize(viewsCount);
        if (storeHeightfield) {
            buffers.resize(viewsCount);
        }
        std::for_each(views.begin(), views.end(), [](const View& view){ assert(view.empty());});

        if (!creatingNew) {
            for (auto iterLat = cellsFrom.latitude; iterLat <= cellsTo.latitude; ++iterLat) {
                for (auto iterLng = cellsFrom.longitude; iterLng <= cellsTo.longitude; ++iterLng) {
                    PlanarIndex iter{iterLat, iterLng};
                    const std::string formattedLatLng = DEMFilenameFormatter::formatLatLng(iter - positivityShift);
                    const std::string filePath
                        = DEMFilenameFormatter::lastVersionNdemPath(normalizedDir(directoryPath), formattedLatLng);
                    MappedFile file{filePath};
                    if (!file.empty()) {
                        ++filesCount;
                        if (verbose) {
                            std::cout << "using " << filePath << std::endl;
                        }
                        auto heightfield = static_cast<Value *>(const_cast<void *>(file.contents()));
                        PlanarIndex viewIndex{iterLat - cellsFrom.latitude, iterLng - cellsFrom.longitude};
                        if (storeHeightfield) {
                            View tmp(heightfield, fileSize.x, fileSize.y);
                            auto index = viewIndexFromCellIndex(viewIndex);
                            views[index] = View(buffers[index], fileSize.x, fileSize.y);
                            tmp.copy(views[index], fileSize.x, fileSize.y);
                        } else {
                            views[viewIndexFromCellIndex(viewIndex)] = View(heightfield, fileSize.x, fileSize.y);
                        }

                        files.push_back(std::move(file)); // TODO: move to else above
                    } else {
                        ++missingFiles;
                        if (verbose) {
                            std::cout << "failed to read " << filePath << " at " << formattedLatLng << std::endl;
                        }
                    }
                }
            }
        }
        if (filesCount == 0 && !storeHeightfield) {
            std::stringstream ss;
            ss << "no files in " << directoryPath << " from " << from << " to " << to;
            throw std::logic_error(ss.str());
        }
        setMargin(0);
        if (verbose && missingFiles > 0) {
            std::cout << "total files: " << filesCount << " missing: " << missingFiles << std::endl;
        }
    }

    void setMargin(Index margin) {
        assert(margin >= 0);
        _indexFrom = PlanarIndex{cellsFrom.latitude * cellsResolution.x + margin, cellsFrom.longitude * cellsResolution.y + margin};
        _indexTo = PlanarIndex{(cellsTo.latitude + 1) * cellsResolution.x - margin, (cellsTo.longitude + 1) * cellsResolution.y - margin};
    }

    void deinit() {
        views.clear();
        files.clear();
    }

    inline bool geodesicIsInside(double latitude, double longitude) const {
        return geodesicIsInside({latitude, longitude});
    }

    inline bool geodesicIsInside(Geodesic geo) const {
        PlanarIndex planar = converter.geodesicToPlanar(geo);
        return planarIsInside(planar);
    }

    inline PlanarIndex geodesicToPlanar(const Geodesic &geo) const {
        return converter.geodesicToPlanar(geo);
    }

    inline Geodesic planarToGeodesic(PlanarIndex planar) const {
        return converter.planarToGeodesic(planar);
    }

    PlanarIndex indexFrom() const {
        return _indexFrom;
    }

    PlanarIndex indexTo() const {
        return _indexTo;
    }

    bool planarIsInside(const PlanarIndex &index) const {
        return index.x < indexTo().x && index.x >= indexFrom().x && index.y < indexTo().y && index.y >= indexFrom().y;
    }

    Vicinity getVicinity(const PlanarIndex &planar) const {
        Vicinity vicinity;
        for (const auto &offsetIndex : offsetsIndexes) {
            auto neighborPlanar = planar + offsetIndex.first;
            if (!planarIsInside(neighborPlanar)) {
                vicinity.edgeMask += (1 << offsetIndex.second);
                continue;
            }
            vicinity.neighbors.emplace_back(SpatialIndex{neighborPlanar, get(neighborPlanar)});
            if (vicinity.neighbors.size() > 1) {
                if (vicinity.neighbors[vicinity.highestNeighbor] < vicinity.neighbors.back()) {
                    vicinity.highestNeighbor = vicinity.neighbors.size() - 1;
                } else if (vicinity.neighbors.back() < vicinity.neighbors[vicinity.lowestNeighbor]) {
                    vicinity.lowestNeighbor = vicinity.neighbors.size() - 1;
                }
            }
        }
        return vicinity;
    }

    uint8_t getEdgeMask(const PlanarIndex &planar) const { // TODO:optimize
        uint8_t mask = 0;
        for (const auto &offsetIndex : offsetsIndexes) {
            if (!planarIsInside(planar + offsetIndex.first)) {
                addToMask(mask, static_cast<uint8_t>(offsetIndex.second));
            }
        }
        return mask;
    }

    inline int viewIndexFromCellIndex(const PlanarIndex &subIndex) const {
        return subIndex.x + subIndex.y * viewsSize.x;
    }

    View& view(PlanarIndex cellIndex) {
        return views[viewIndexFromCellIndex(cellIndex)];
    }

    static constexpr Value kFallBackValue() {
        return kOceanFloor;
    }

    Value get(const PlanarIndex &index, bool dbg = false) const {
        if (!planarIsInside(index)) {
            return kFallBackValue();
        }
        PlanarIndex subIndex = PlanarIndex{index.x % cellsResolution.x, index.y % cellsResolution.y};
        uint16_t viewIndex = viewIndexFromCellIndex({index.x / cellsResolution.x, index.y / cellsResolution.y});

        auto &view = views.at(viewIndex);

        if ((index - debugPlanar).manhattanLength() < 2) {
            std::cout << "dbg get " << index << " vi: " << viewIndex << " subIndex: " << subIndex << " view size: " << view.size << " view.cosize: " << view.cosize << endl;
        }

        if (view.empty()) {
            return kFallBackValue();
        } else {
            return view.get(subIndex);
        }
    }

    Value& getRef(const PlanarIndex &index) {
        assert(isStoringField());
        PlanarIndex subIndex = views.size() == 1 ? index : PlanarIndex{index.x % cellsResolution.x, index.y % cellsResolution.y};
        uint16_t viewIndex = views.size() == 1 ? 0 : viewIndexFromCellIndex({index.x / cellsResolution.x, index.y / cellsResolution.y});
        auto &view = views.at(viewIndex);

        if (view.empty()) {
            view = View(buffers[viewIndex], fileSize.x, fileSize.y);
        }

        return view.get(subIndex);
    }

    void set(const PlanarIndex &index, Value value) {
        if (!isStoringField()) {
            throw std::logic_error("trying to set value in non field storing grid");
        }
        PlanarIndex subIndex{index.x % cellsResolution.x, index.y % cellsResolution.y};
        uint16_t viewIndex = viewIndexFromCellIndex({index.x / cellsResolution.x, index.y / cellsResolution.y});

        auto &view = views.at(viewIndex);
        if (view.empty()) {
//            std::cout << "creating view at " << viewIndex << endl;
            view = View(buffers.at(viewIndex), fileSize.x, fileSize.y);
        }
        view.get(subIndex) = value;
    }

    bool isStoringField() const {
        return !buffers.empty();
    }

    bool empty() const {
        for (auto &view : views) {
            if (!view.empty()) {
                return false;
            }
        }
        return true;
    }

    void doForEachPlanar(const std::function<void(const PlanarIndex& planar)> &callback) const {
        for (Index la = indexFrom().x; la < indexTo().x; ++la) {
            for (Index lo = indexFrom().y; lo < indexTo().y; ++lo) {
                callback(PlanarIndex{la, lo});
            }
        }
    }

    size_t planarIndexSize() const {
        PlanarIndex size = indexTo() - indexFrom();
        return size.x * size.y;
    }

    void serialize(const std::string &outputDir, const std::string &extension) const {
        assert(isStoringField());
        for (Index la = cellsFrom.latitude; la <= cellsTo.latitude; ++la) {
            for (Index lo = cellsFrom.longitude; lo <= cellsTo.longitude; ++lo) {
                PlanarIndex cellIndex{la , lo};
                uint16_t viewIndex = viewIndexFromCellIndex(cellIndex);
                if (views[viewIndex].empty()) {
                    continue;
                }
                PlanarIndex cellOrigin = cellIndex - positivityShift;
                auto path = outputDir + DEMFilenameFormatter::formatLatLng(cellOrigin) + "." + extension;
                assert(!buffers[viewIndex].empty());
                serializeNonEmptyVector(buffers[viewIndex], path);
            }
        }
    }
};

template <typename View = FieldView>
struct DEMGrid : PlanarGrid<View, Altitude> {
    using BaseGrid = typename DEMGrid::PlanarGrid;
    using BaseGrid::planarIsInside;
    using BaseGrid::get;
    using BaseGrid::BaseGrid;

    mutable struct {
        SpatialIndex highestNeighbor;
        PlanarIndexDir highestNeighborDir;
        uint8_t dirsUpMask; // TODO: union to speed up?
        uint8_t edgeMask;
    } runningVicinityInfo;

    mutable SpatialIndex runningHighestNeighbor;

    bool isPeak(const PlanarIndex &planar, const Altitude &altitude, bool strict = false, bool flipTheWorld = false) const {
        unsigned aboveCount = 0;
        int dir = flipTheWorld ? -1 : 1;
        for (const auto &neighbor : immediateVicinity(planar)) {
            if (!planarIsInside(neighbor)) {
                ++aboveCount;
                continue;
            }
            auto neighborAltitude = BaseGrid::get(neighbor);
            if (neighborAltitude * dir > altitude * dir) {
                return false;
            } else if (neighborAltitude * dir < altitude * dir) {
                ++aboveCount;
            }
        }
        return strict ? aboveCount == kMaxNeighborsCount : aboveCount > 1;
    }

    SpatialIndex getHighestNeighbour(SpatialIndex &spatial) const {
        auto highestNeighbour = spatial;
        for (const auto &planar : immediateVicinity(spatial.planar)) {
            SpatialIndex neighbour{planar, get(planar)};
            if (highestNeighbour < neighbour) {
                std::swap(highestNeighbour, neighbour);
            }
        }
        return highestNeighbour;
    }

    bool isPeakDefinitevly(const PlanarIndex &planar) const {
        SpatialIndex spatial{planar, get(planar)};
        return isPeakDefinitevly(spatial);
    }

    bool isPeakDefinitevly(const SpatialIndex &spatial) const {
        for (const auto &neighbor : immediateVicinity(spatial.planar)) {
            if (!planarIsInside(neighbor)) {
                continue;
            }
            SpatialIndex neighborSpatial{neighbor, get(neighbor)};
            if (spatial < neighborSpatial) {
                return false;
            }
        }
        return true;
    }

    bool definitivelyLower(const PlanarIndex &first, const PlanarIndex &second) const {
        Altitude firstAlt = get(first);
        Altitude secondAlt = get(second);
        if (firstAlt == secondAlt) {
            return first < second;
        } else {
            return firstAlt < secondAlt;
        }
    }

    template <uint16_t radius>
    bool isPotentialCol(const SpatialIndex &spatial, bool updateRunningVicinityInfo = false) const {
        bool verbose = false;
        static std::vector<PlanarIndex> potentialRidgeStarts;
        potentialRidgeStarts.clear();
        if (verbose){
            std::cout << "checking for col: " << spatial << endl;
        }
        if (updateRunningVicinityInfo) {
            runningVicinityInfo.dirsUpMask = 0;
            runningVicinityInfo.edgeMask = 0;
            runningVicinityInfo.highestNeighborDir = kPlanarIndexNoDir;
        }

        for (uint8_t offsetIndex = 0; offsetIndex < kMaxNeighborsCount; ++offsetIndex) {
            auto offsetDir = PlanarIndexDir(offsetIndex);
            PlanarIndex neighborPlanar = spatial.planar + kNeighborOffsets[offsetIndex];

            if (updateRunningVicinityInfo) {
                if (!planarIsInside(neighborPlanar)) {
                    addToMask(runningVicinityInfo.edgeMask, offsetIndex);
                    assert(checkMask(runningVicinityInfo.edgeMask, offsetIndex));
                } else {
                    assert(!checkMask(runningVicinityInfo.edgeMask, offsetIndex));
                }
            }
            auto nspatial = SpatialIndex{neighborPlanar, get(neighborPlanar)};

            if (spatial.planar == debugPlanar) {
                std::cout << "dir: " << offsetDir << " spatial: " << nspatial << endl;
            }

            if (verbose) {
                std::cout << "n: " << nspatial << endl;
            }

            if (spatial < nspatial) {
                potentialRidgeStarts.push_back(neighborPlanar);
                if (updateRunningVicinityInfo) {
                    if (spatial.planar == debugPlanar) {
                        std::cout << "adding to dirsUp\n";
                    }
                    addToMask(runningVicinityInfo.dirsUpMask, offsetIndex);
                    if (runningVicinityInfo.highestNeighborDir == kPlanarIndexNoDir || runningVicinityInfo.highestNeighbor < nspatial) {
                        runningVicinityInfo.highestNeighborDir = offsetDir;
                        runningVicinityInfo.highestNeighbor = nspatial;
                    }
                }
            }
        }
        if (potentialRidgeStarts.size() < 2) {
            if (verbose) {
                std::cout << "not a col\n";
            }
            return false;
        }

        static VicinityBitfield<radius> firstRidgeComponent;
        firstRidgeComponent.clear();
        static VicinityBitfield<radius> allProcessedPlanars;
        allProcessedPlanars.clear();
        allProcessedPlanars.set({0, 0});
        static std::vector<PlanarIndex> frontier;
        frontier.clear();
        frontier.push_back(potentialRidgeStarts.front());
        allProcessedPlanars.set(spatial.planar - potentialRidgeStarts.front());
        uint16_t immediateNeighborsInComponent = 1;
        while (!frontier.empty()) {
            auto next = frontier.back();
            if (verbose) {
                std::cout << "next: " << next << endl;
            }
            auto nextOffset = spatial.planar - next;
            firstRidgeComponent.set(nextOffset);

            frontier.pop_back();
            for (const auto &neighbor : immediateVicinity(next)) {
                if (!planarIsInside(neighbor)) {
                    continue;
                }
                PlanarIndex neighborOffset = spatial.planar - neighbor;
                if (!allProcessedPlanars.isInside(neighborOffset)) {
                    if (verbose) {
                        std::cout << "not in vicinity: " << neighborOffset << endl;
                    }
                    continue;
                }
                if (allProcessedPlanars.get(neighborOffset)) {
                    if (verbose) {
                        std::cout << "already processed: " << neighborOffset << std::endl;
                    }
                    continue;
                }

                SpatialIndex neighborSpatial{neighbor, get(neighbor)};
                if (spatial < neighborSpatial) {
                    if (verbose) {
                        std::cout << "adding to frontier: " << neighborOffset << endl;
                    }
                    if (neighborOffset.manhattanLength() == 1) {
                        ++immediateNeighborsInComponent;
                        if (immediateNeighborsInComponent == potentialRidgeStarts.size()) {
                            return false;
                        }
                    }
                    frontier.push_back(neighbor);
                    allProcessedPlanars.set(neighborOffset);
                } else {
                    if (verbose) {
                        std::cout << "not adding to frontier: " << neighborOffset << " = " << neighborSpatial << endl;
                    }
                    allProcessedPlanars.set(neighborOffset);
                }
            }

        }

        if (spatial.planar == debugPlanar) {
            std::cout << "debug planar is a col\n";
        }

        return true;
    }

};

template <typename View, uint16_t kLogValueBits>
struct CompactGrid : PlanarGrid<View, uint16_t> {
    using ContainerValue = uint16_t;
    static constexpr uint16_t kContainerLogBits = 4;
    static constexpr uint16_t kContainerBits = 1 << kContainerLogBits;
    static constexpr uint16_t kValueBits = 1 << kLogValueBits;
    static constexpr uint16_t kValueMax = 1 << kValueBits;
    static constexpr uint16_t kLogValuesInContainer = kContainerLogBits - kLogValueBits;
    static constexpr uint16_t kValuesInContainer = 1 << kLogValuesInContainer;
    static constexpr ContainerValue kSubMask = (1 << kLogValuesInContainer) - 1;
    static constexpr uint16_t kValueMask = (1 << kValueBits) - 1;
    using BaseGrid = typename CompactGrid::PlanarGrid;

    CompactGrid(const std::string &directoryPath,
                const PlanarIndex &from,
                const PlanarIndex &to,
                const PlanarIndex &fileSize,
                PlanarIndex cellsResolution,
                bool verbose = false,
                bool storeHeightfield = false)
    : BaseGrid(directoryPath, from, to, {fileSize.x >> kLogValuesInContainer, fileSize.y}, {cellsResolution.x >> kLogValuesInContainer, cellsResolution.y}, verbose, storeHeightfield) {
        static_assert(sizeof(ContainerValue) * CHAR_BIT == kContainerBits, "inconsistent bits for container value in compact grid");
        assert(kValueMask == 3);
        assert(kSubMask == 7);
        std::cout << "Compact grid: " << sizeof(ContainerValue) << " " << kValueBits << " " << kLogValuesInContainer << endl;
    }

    uint16_t get(const PlanarIndex &planar, bool dbg = false) const {
        PlanarIndex basePlanar = getBasePlanar(planar);

        auto containerValue = BaseGrid::get(basePlanar);
        auto subIndex = planar.x & kSubMask;
        if (dbg) {
            BaseGrid::get(basePlanar, true);
            std::cout << "get cv = " << containerValue << " si = " << subIndex << " basePlanar: " << basePlanar << endl;
        }

        assert(subIndex < kValuesInContainer);
        return (containerValue >> (subIndex * kValueBits)) & kValueMask;
    }

    void set(const PlanarIndex &planar, uint16_t value) {
        assert(value < kValueMax);
//        auto dbgValue = get(debugPlanar);
        PlanarIndex basePlanar = getBasePlanar(planar);
        auto containerValue = BaseGrid::get(basePlanar); // TODO: getRef
//        auto containerValueBefore = containerValue;
        auto subIndex = planar.x & kSubMask;
        assert(subIndex < kValuesInContainer);
        ContainerValue cleanMask = ~(kValueMask << (subIndex * kValueBits));
        containerValue &= cleanMask;
        auto setMask = value << (subIndex * kValueBits) ;
        containerValue |= setMask;

        BaseGrid::set(basePlanar, containerValue);

    }
    
    static PlanarIndex getBasePlanar(const PlanarIndex &planar) {
        return {planar.x >> kLogValuesInContainer, planar.y};
    }
    
    uint8_t getEdgeMask(const PlanarIndex &planar) const {
        return BaseGrid::getEdgeMask(getBasePlanar(planar));
    }
    
    uint8_t planarIsInside(const PlanarIndex &planar) const {
        return BaseGrid::planarIsInside(getBasePlanar(planar));
    }
    
};

using RegularDEMGrid = DEMGrid<FieldView>;
using TiledDEMGrid = DEMGrid<TiledFieldView>;

struct SDFileFormat {
    using DEM = RegularDEMGrid;
    static constexpr PlanarIndex fileSize{1200, 1200};
    static constexpr PlanarIndex cellsResolution{1200, 1200};
    static std::string prefix() {
        return "sd_";
    }
    static DEM getDEM(const std::string &demsDirectory, const Rectangle &region) {
        return DEM{demsDirectory, region.origin, region.end - kDiagonalShift, {fileSize}, cellsResolution};
    }
};

struct HDFileFormat {
    using DEM = TiledDEMGrid;
    static constexpr PlanarIndex fileSize{3616, 3648};
    static constexpr PlanarIndex cellsResolution{3600, 3600};
    static std::string prefix() {
        return "hd_";
    }

    static DEM getDEM(const std::string &demsDirectory, const Rectangle &region) {
        return DEM{demsDirectory, region.origin, region.end - kDiagonalShift, {fileSize}, cellsResolution};
    }
};


//TiledDEMGrid&& buildTiledGrid(const std::string &demsDirectory, const PlanarIndex &from, const PlanarIndex &to, bool verbose = false) {
//    return std::move(TiledDEMGrid(demsDirectory, from, to, {kBigTiledFileSizeX, kBigTiledFileSizeY}, 3600, verbose));
//}
