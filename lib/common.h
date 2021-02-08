#pragma once
#undef NDEBUG

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>
#include <filesystem>

#ifdef DEBUG_ASSERTS
# define debug_assert(EX) assert(EX)
#else
# define debug_assert(EX)
#endif

#ifdef QTREE_DEBUG
#define qtree_debug(NX) NX
#else
#define qtree_debug(NX)
#endif

using Path = std::filesystem::path;

template <typename T>
inline std::string outstrImpl(const std::string &name, const T &t) {
    return name + ": " + std::to_string(t);
}

inline std::string outstrImpl(const std::string &name, const std::string &t) {
    return name + ": " + t;
}

#define OUTSTR(x) (outstrImpl(#x, (x)) + " ")

using Index = int32_t;
using Altitude = int16_t;
using Distance = double;
using IslandIndex = int32_t;
using Degrees = double;
using PeakId = uint64_t;
using namespace std::literals;
constexpr auto endl = "\n";

static constexpr Altitude kNoAltitude = -1235;
static const Altitude kOceanFloor = 0;
static const Degrees kNoDegrees = 9999;
static const size_t kNoIndexSizeT = std::numeric_limits<size_t>::max();
static const uint32_t k32NoIndex = std::numeric_limits<uint32_t>::max();
static const Index kNoIndex = std::numeric_limits<Index>::max();

static const IslandIndex kNoIsland = -1;
static const int kBigTiledFileSizeX = 3616;
static const int kBigTiledFileSizeY = 3648;
static const int kHDDemSize = 3600;

static const Altitude kHighNansAltitude = 8851;
static const Altitude kLowNansAltitude = -500;
constexpr PeakId kNoId = 0;

struct Geodesic {
    union {
        struct {
            double latitude, longitude;
        };
        struct {
            double x, y;
        };
        struct {
            double c[2];
        };
    };

    bool isNormalized() const {
        return std::abs(latitude) <= 90 && std::abs(longitude) <= 180;
    }

    void normalize() {
        if (isNormalized()) {
            return;
        }
        while (std::abs(latitude) > 180) {
            latitude += latitude > 0 ? -360 : 360;
        }
        if (std::abs(latitude) > 90) {
            latitude = (latitude > 0 ? 1 : -1) * 180 - latitude;
            longitude += 180;
        }
        while (longitude > 180 || longitude <= -180) {
            longitude += longitude > 0 ? -360 : 360;
        }
        if (longitude == -180.) {
            longitude = 180;
        }
    }

    inline bool operator <(const Geodesic &other) const {
        return x < other.x || (x == other.x && y < other.y);
    }

    inline bool operator==(const Geodesic &other) const {
        return x == other.x && y == other.y;
    }

    std::string str() const {
        return "(" + std::to_string(latitude) + ", " + std::to_string(longitude) + ")";
    }
};

inline std::string normalizedDir(const std::string &dir) {
    return dir.back() == '/' ? dir : dir + '/';
}
