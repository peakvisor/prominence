#pragma once

#include "common.h"

#include <iostream>
#include <array>
#include <functional>
#include <map>
#include <unordered_set>

struct Rectangle;

struct PlanarIndex {
    static PlanarIndex outputZero;
    
    union {
        struct {
            Index x, y;
        };
        struct {
            Index latitude, longitude; // 'vertical' latitude = 'x'
        };
        struct {
            Index c[2];
        };
        long prehash;
    };

    template <int sizeLog2>
    constexpr Index linearize_() const {
        return x + (y << sizeLog2);
    }
    
    constexpr Index linearize(Index size) const {
        return x + y * size;
    }

    constexpr int64_t linearizeBig(Index size) const {
        return x + y * static_cast<int64_t>(size);
    }

    constexpr PlanarIndex &operator +=(const PlanarIndex &other) {
        x += other.x; y += other.y;
        return *this;
    }
    
    constexpr PlanarIndex operator +(const PlanarIndex &other) const {
        return {x + other.x, y + other.y};
    }
    
    constexpr PlanarIndex operator -(const PlanarIndex &other) const {
        return {x - other.x, y - other.y};
    }
    
    constexpr PlanarIndex operator -() const {
        return {-x, -y};
    }
    
    inline PlanarIndex operator *(Index scale) const {
        return {x * scale, y * scale};
    }
    
    inline PlanarIndex operator /(Index scale) const {
        return {x / scale, y / scale};
    }
    
    inline bool operator <(const PlanarIndex &other) const {
        return x == other.x ? y < other.y : x < other.x;
    }
    
    constexpr bool operator ==(const PlanarIndex &other) const {
        return x == other.x && y == other.y;
    }
    
    constexpr bool operator !=(const PlanarIndex &other) const {
        return x != other.x || y != other.y;
    }
    
    inline static Index nonnegativeResidual(Index index, Index modulus) {
        auto residual = index % modulus;
        residual += residual < 0 ? modulus : 0;
        return residual;
    }
    
    [[nodiscard]] inline PlanarIndex nonnegativeResidual(Index modulus) const {
        return {nonnegativeResidual(x, modulus), nonnegativeResidual(y, modulus)};
    }
    
    [[nodiscard]] inline constexpr uint64_t zHashLike() const {
        return ((static_cast<uint64_t>(x) >> 4u) << 32u) | ((static_cast<uint64_t>(y) >> 4u));
    }
    
    [[nodiscard]] inline bool zHashLikeCompare(const PlanarIndex &other) const {
        return zHashLike() < other.zHashLike();
    }
    
    [[nodiscard]] inline bool isInsideRectangleWithSize(
            const PlanarIndex &origin, const PlanarIndex &size) const
    {
        PlanarIndex diff = {x - origin.x, y - origin.y};
        return 0 <= diff.x && diff.x < size.x && 0 <= diff.y && diff.y < size.y;
    }

    [[nodiscard]] inline bool isInsideRectangle(const PlanarIndex &origin, const PlanarIndex &end) const {
        return origin.x <= x && x < end.x && origin.y <= y && y < end.y;
    }
    
    [[nodiscard]] bool isInsideRectangle(const Rectangle &rectangle) const;
    
    [[nodiscard]] inline bool isInsideSquare(const PlanarIndex &origin, const Index &side) const {
        PlanarIndex diff = {x - origin.x, y - origin.y};
        return 0 <= diff.x && diff.x < side && 0 <= diff.y && diff.y < side;
    }
    
    [[nodiscard]] inline Index manhattanLength() const {
        return std::abs(x) + std::abs(y);
    }
};

inline PlanarIndex planarFromIndexedCoordinates(int firstCoordinateIndex,
        Index firstCoordinate, Index secondCoordinate)
{
    PlanarIndex res;
    res.c[firstCoordinateIndex] = firstCoordinate;
    res.c[1 - firstCoordinateIndex] = secondCoordinate;
    return res;
}

inline PlanarIndex planarFromStr(const std::string &laStr, const std::string &loStr) {
    return {std::stoi(laStr), std::stoi(loStr)};
}

struct PlanarIndexHash {
    inline size_t operator ()(const PlanarIndex &planar) const {
        return planar.prehash;
    }
};

using PlanarSet = std::unordered_set<PlanarIndex, PlanarIndexHash>;

struct Rectangle {
    PlanarIndex origin;
    PlanarIndex end;
    
    constexpr uint32_t area() const {
        auto res = (end.x - origin.x) * (end.y - origin.y);
        assert(res >= 0);
        return res;
    }
    
    constexpr Rectangle inflated(Index inflationSize) const {
        return Rectangle{origin - PlanarIndex{inflationSize, inflationSize}, end + PlanarIndex{inflationSize, inflationSize}};
    }
    
    constexpr void boundBy(const Rectangle &boundingRectangle) {
        for (int i = 0; i < 2; ++i) {
            if (origin.c[i] < boundingRectangle.origin.c[i]) {
                origin.c[i] = boundingRectangle.origin.c[i];
            }
            if (end.c[i] > boundingRectangle.end.c[i]) {
                end.c[i] = boundingRectangle.end.c[i];
            }
        }
    }
    
    void applyToAll(std::function<void(const PlanarIndex&)> &lambda) const {
        for (Index la = origin.x; la < end.x; ++la) {
            for (Index lo = origin.y; lo < end.y; ++lo) {
                lambda({la, lo});
            }
        }
    }
};

constexpr PlanarIndex min(const PlanarIndex &lhs, const PlanarIndex &rhs) {
    return {std::min(lhs.x, rhs.x), std::min(lhs.y, rhs.y)};
}

constexpr PlanarIndex max(const PlanarIndex &lhs, const PlanarIndex &rhs) {
    return {std::max(lhs.x, rhs.x), std::max(lhs.y, rhs.y)};
}

constexpr bool rectanglesIntersect(const Rectangle &first, const Rectangle &second) {
    for (int i = 0; i < 2; ++i) {
        if (first.end.c[i] <= second.origin.c[i] || first.origin.c[i] >= second.end.c[i]) {
            return false;
        }
    }
    return true;
}

constexpr Rectangle rectangularShell(const Rectangle &lhs, const Rectangle &rhs) {
    return {min(lhs.origin, rhs.origin), max(lhs.end, rhs.end)};
}

constexpr bool rectanglesNear(const Rectangle &first, const Rectangle &second) {
    return rectanglesIntersect(first.inflated(1), second);
}

static constexpr PlanarIndex kLatLngFrom{-90, -180},
                             kLatLngTo{90, 180};


struct SpatialIndex {
    PlanarIndex planar;
    Altitude altitude;
    
    inline bool operator<(const SpatialIndex& other) const {
        return altitude == other.altitude ? planar < other.planar : altitude < other.altitude;
    }
        
    inline bool operator>(const SpatialIndex& other) const {
        return other < *this;
    }
    
    inline bool operator==(const SpatialIndex& other) const {
        return planar == other.planar && altitude == other.altitude;
    }

    inline bool operator!=(const SpatialIndex& other) const {
        return planar != other.planar || altitude != other.altitude;
    }
};

struct LinearIndexer {
    inline int operator()(const PlanarIndex &index, int size) const {
        return index.linearize(size);
    }
    inline PlanarIndex reverse(int linear, int size) {
        return {linear % size, linear / size};
    }
};

template <uint16_t tileLog2>
struct TiledIndexer {
    static constexpr int twiceTileLog2() { return tileLog2 << 1u; }
    static constexpr int tileMask() { return (1u << tileLog2) - 1; }
    inline int operator()(const PlanarIndex &index, int size) const {
        auto superIndex = PlanarIndex{index.x >> tileLog2, index.y >> tileLog2}.linearize(size >> tileLog2);
        return PlanarIndex{index.x & tileMask(), index.y & tileMask()}.linearize_<tileLog2>() | (superIndex << twiceTileLog2());
    }
};

template <uint16_t xTileLog2, uint16_t yTileLog2>
struct AssymetricTiledIndexer {
    static constexpr uint16_t superIndexResolutionLog2() { return xTileLog2 + yTileLog2; }
    static constexpr int xTileMask() { return (1u << xTileLog2) - 1; }
    static constexpr int yTileMask() { return (1u << yTileLog2) - 1; }
    inline int operator()(const PlanarIndex &index, int size) const {
        auto superIndex = PlanarIndex{index.x >> xTileLog2, index.y >> yTileLog2}.linearize(size >> xTileLog2);
        return PlanarIndex{index.x & xTileMask(), index.y & yTileMask()}.linearize_<xTileLog2>() | (superIndex << superIndexResolutionLog2());
    }
    inline PlanarIndex reverse(int linear, int size) {
        int tileArea = 1u << superIndexResolutionLog2();
        PlanarIndex tileSize = {1u << xTileLog2, 1u << yTileLog2};
        int superSize = size >> xTileLog2;
        int superLinear = linear / tileArea;
        PlanarIndex superPlanar = {superLinear % superSize, superLinear / superSize};
        PlanarIndex superOrigin = {superPlanar.x * tileSize.x, superPlanar.y * tileSize.y};
        
        int subLinear = linear % tileArea;
        PlanarIndex subPlanar = {subLinear % (1u << xTileLog2), subLinear / (1u << xTileLog2)};
        return superOrigin + subPlanar;
    }
};

static constexpr PlanarIndex kDiagonalShift{1, 1};

//std::ostream& operator<<(std::ostream& out, const PlanarIndex& index);
inline std::ostream& operator<< (std::ostream& out, const PlanarIndex& index) {
    out << "(" << index.x - PlanarIndex::outputZero.x << ", " << index.y  - PlanarIndex::outputZero.y << ")";
    return out;
}

inline std::ostream& operator<<(std::ostream& out, const SpatialIndex& index) {
    out << "(" << index.planar.x << ", " << index.planar.y << ", " << index.altitude << ")";
    return out;
}

inline std::ostream& operator<< (std::ostream& out, const Geodesic &geo) {
    out << "(" << geo.latitude << ", " << geo.longitude << ")";
    return out;
}

// we are using: k(Full)NeighborOffsets[i ^ 1] == -k(Full)NeighborOffsets[i]
static constexpr PlanarIndex kNeighborOffsets[] = {
    {1, 0}, {-1, 0}, {0, 1}, {0, -1}
};

enum PlanarIndexDir : uint8_t {
    kPlanarIndexRight = 0,
    kPlanarIndexLeft,
    kPlanarIndexUp,
    kPlanarIndexDown,
    kPlanarIndexNoDir
};

constexpr PlanarIndex kFullNeighborsOffsets[] = {
    {1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {-1, -1}, {1, -1}, {-1, 1}
};

constexpr PlanarIndex kFullNeighborhoodTraversal[] = {
    {1, 0}, {1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}
};

constexpr unsigned kMaxNeighborsCount = sizeof(kNeighborOffsets) / sizeof(kNeighborOffsets[0]);
constexpr unsigned kMaxFullNeighborsCount = sizeof(kFullNeighborsOffsets) / sizeof(kFullNeighborsOffsets[0]);

inline std::array<PlanarIndex, kMaxNeighborsCount> immediateVicinity(const PlanarIndex &planar) {
    return {planar + kNeighborOffsets[0], planar + kNeighborOffsets[1],
        planar + kNeighborOffsets[2], planar + kNeighborOffsets[3]};
}
std::array<PlanarIndex, kMaxNeighborsCount> neighborsScaledOffsets(int scale);

std::array<PlanarIndex, kMaxFullNeighborsCount> fullImmediateVicinityTraversal(const PlanarIndex &planar);

static const PlanarIndex debugPlanar{1733999, 24485};

static const std::map<PlanarIndex, PlanarIndexDir> offsetsIndexes = // TODO: unordered?
    {{kNeighborOffsets[kPlanarIndexRight], kPlanarIndexRight},
     {kNeighborOffsets[kPlanarIndexLeft], kPlanarIndexLeft},
     {kNeighborOffsets[kPlanarIndexUp], kPlanarIndexUp},
     {kNeighborOffsets[kPlanarIndexDown], kPlanarIndexDown}};

uint8_t offsetMask(PlanarIndexDir offsetDir);
PlanarIndexDir offsetDirFromMask(uint8_t mask);
void addToMask(uint8_t &mask, uint8_t offsetIndex);
bool checkMask(uint8_t mask, uint8_t offsetIndex);
PlanarIndex nextInDir(const PlanarIndex &planar, PlanarIndexDir dir);
