#pragma once
#include "common.h"

#include <cmath>
#include <utility>
#include <array>
#include <sstream>
#include <iomanip>
#include <algorithm>

constexpr double kEarthRadius = 6378.137;
constexpr double kPiDegrees = 180.;
constexpr double kOneDegreeInRads = M_PI / 180.;
constexpr double kOneRadInDegrees = 1 / kOneDegreeInRads;

struct MinMax {
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
};

struct Cartesian {
    double x, y, z;

    Cartesian operator*(double a) const {
        return Cartesian{x * a, y * a, z * a};
    }

    Cartesian operator+(const Cartesian &o) const {
        return {x + o.x, y + o.y, z + o.z};
    }

    Cartesian operator-(const Cartesian &o) const {
        return {x - o.x, y - o.y, z - o.z};
    }

    double operator*(const Cartesian &o) const {
        return x * o.x + y * o.y + z * o.z;
    }

    double length() const {
        return sqrt(x * x + y * y + z * z);
    }

    void normalize() {
        auto l = length();
        x = x / l;
        y = y / l;
        z = z / l;
    }
};

inline Cartesian geoToCartesian(double latitude, double longitude) {
    return {std::cos(latitude) * std::cos(longitude), std::cos(latitude) * std::sin(longitude), std::sin(latitude)};
}

inline Cartesian meridianNormal(double longitude) {
    return {std::sin(longitude), -std::cos(longitude), 0.};
}

inline double longitudeFromCartesian(const Cartesian &v) {
    return std::atan(v.y / v.x) - M_PI * (v.y < 0 && v.x < 0) + M_PI * (v.x < 0 && v.y >= 0);
}

inline double latitudeFromCartesian(const Cartesian &v) {
    return std::asin(v.z);
}

constexpr double degreesFromRads(double rads) {
    return rads * kOneRadInDegrees;
}

constexpr double radsFromDegrees(double degrees) {
    return degrees * kOneDegreeInRads;
}

constexpr Geodesic radsFromDegrees(const Geodesic &degrees) {
    return {radsFromDegrees(degrees.latitude), radsFromDegrees(degrees.longitude)};
}

inline std::ostream& operator<<(std::ostream &out, const Cartesian &c) {
    out << "(" << c.x << ", " << c.y << ", " << c.z << ")";
    return out;
}

constexpr double clamp(double value, double maxAbs = 1) { // too fast to check that maxAbs > 0
    return std::min(std::max(value, -maxAbs), +maxAbs);
}

struct DistanceTools {
    static constexpr double kHalfPi = M_PI_2;
    static constexpr double kProjEps = 1.0e-12;
    static constexpr double kSwitchedLongitudeEps = 2.0e-6;
    static constexpr double kRadEps = 1.0e-12;
    static constexpr double kMaxDegreesLatitude = kPiDegrees / 2;
    static constexpr double kMaxDegreesLongitude = kPiDegrees;
    static constexpr double kMaxRads = M_PI;
    static constexpr double kMinRads = -kMaxRads;

    static void normalizeRad(double &rad) {
        while (rad > kMaxRads) {
            rad -= 2 * M_PI;
        }
        while (rad < kMinRads) {
            rad += 2 * M_PI;
        }
    }

    static double addRads(double x, double y) {
        double res = x + y;
        normalizeRad(res);
        return res;
    }

    static Cartesian closestCartesianOnMeridian(const Cartesian& from, double meridian) {

        auto normal = meridianNormal(meridian);

        double fromNormalProduct = from * normal;
        if (std::abs(1 - std::abs(fromNormalProduct)) < kProjEps) {
            return {std::cos(meridian), std::sin(meridian), 0};
        }
        auto result = (from - normal * fromNormalProduct);
        result.normalize();
        return result;
    }

    static double minLatDiff(double from, double toBottom, double toTop) {
        if (from > toTop) {
            return from - toTop;
        } else if (from < toBottom) {
            return toBottom - from;
        } else {
            return 0.;
        }
    }

    static MinMax minMaxCircleDiff(double from,
                            double toLeft, double toRight,
                            double &argmin, double &argmax) {
        normalizeRad(from);
        normalizeRad(toLeft);
        normalizeRad(toRight);
        if (toLeft > toRight) {
            std::swap(toLeft, toRight);
        }
        MinMax res;
        if (from < toLeft || from > toRight) {
            auto distanceToLeft = toLeft - from;
            if (distanceToLeft < 0) {
                distanceToLeft += 2 * M_PI;
            }
            auto distanceToRight = from - toRight;
            if (distanceToRight < 0) {
                distanceToRight += 2 * M_PI;
            }

            argmin = distanceToRight < distanceToLeft ? toRight : toLeft;
            res.min = std::min(distanceToRight, distanceToLeft);

            res.max = res.min + toRight - toLeft;
            argmax = distanceToRight > distanceToLeft ? toRight : toLeft;
        } else {
            argmin = from;
            res.min = 0.;

            auto distanceToLeft = from - toLeft;
            auto distanceToRight = toRight - from;

            argmax = distanceToRight > distanceToLeft ? toRight : toLeft;
            res.max = std::max(distanceToRight, distanceToLeft);
        }

        if (res.max > M_PI) {
            res.max = M_PI;
            argmax = addRads(from, M_PI);
            if (argmax > toRight + kRadEps || argmax < toLeft - kRadEps) {
                std::cout << "from: " << degreesFromRads(from)
                          << " to [" << degreesFromRads(toLeft) << ", " << degreesFromRads(toRight)
                          << "], argmax = \n" << degreesFromRads(argmax) << std::endl;
                assert(false);
            }
        }
        return res;
    }

    static MinMax minMaxCircleDiff(double from,
                            double toLeft, double toRight) {
        double argmin, argmax;
        return minMaxCircleDiff(from, toLeft, toRight, argmin, argmax);
    }

    static double maxLatDiff(double from, double toBottom, double toTop) {
        return std::max(std::abs(from - toBottom), std::abs(from - toTop));
    }

    static double distanceToSegmentOnMeridian(double fromLat, double fromLong, double meridian,
                                              double toLatBottom, double toLatTop, bool min) {
        Cartesian from = geoToCartesian(fromLat, fromLong);
        bool noProjectionNeeded = std::abs(meridian - fromLong) < kProjEps || std::abs(std::abs(meridian - fromLong) - M_PI) < kProjEps;
        Cartesian proj = noProjectionNeeded ? from : closestCartesianOnMeridian(from, meridian);
        double projLat = latitudeFromCartesian(proj);
        double projLong = std::abs(std::abs(projLat) - kHalfPi) < kProjEps ? meridian : longitudeFromCartesian(proj);
        double longDiff = std::abs(projLong - meridian);

        double diffFromPi = std::abs(longDiff - M_PI);
        bool switchedLongitude = diffFromPi < kSwitchedLongitudeEps;
        if (!(switchedLongitude || longDiff < kSwitchedLongitudeEps || std::abs(longDiff - 2 * M_PI) < kSwitchedLongitudeEps)) {
            std::stringstream longDiffStream;
            longDiffStream << std::setprecision(20) << longDiff << " " << std::setprecision(20) << longDiff - M_PI;
            longDiffStream << "\nargs: " << fromLat << ", " << fromLong << ", " << meridian << ", "
                      << toLatBottom << ", " << toLatTop << ", " << min << std::endl;
            throw std::logic_error("inconsistent switched longitude longDiff " + longDiffStream.str());
        };

        double effectiveProjLat = switchedLongitude ? M_PI - projLat : projLat;
        normalizeRad(effectiveProjLat);
        double argmin, argmax;

        auto minMax = minMaxCircleDiff(effectiveProjLat, toLatBottom, toLatTop, argmin, argmax);
        double latDiff = min ? minMax.min : minMax.max;
        double firstCos = noProjectionNeeded ? 1. : from * proj;
        double targetCos = firstCos * std::cos(latDiff);
        double distance;
        if (targetCos >= 1) {
            distance = 0;
        } else if (targetCos <= -1) {
            distance = M_PI;
        } else {
            distance = std::acos(targetCos);
        }

        return distance;
    }

    static MinMax minMaxDistanceToQuadrant(double fromLat, double fromLong,
                                    double toLatBottom, double toLatTop,
                                    double toLongLeft, double toLongRight) {
        assert(toLatTop >= toLatBottom);
        assert(toLongRight >= toLongLeft);

        double closestLong;
        double farthestLong;

        minMaxCircleDiff(fromLong, toLongLeft, toLongRight, closestLong, farthestLong);

        return {distanceToSegmentOnMeridian(fromLat, fromLong, closestLong, toLatBottom, toLatTop, true) * kEarthRadius,
                distanceToSegmentOnMeridian(fromLat, fromLong, farthestLong, toLatBottom, toLatTop, false) * kEarthRadius};
    }

    static MinMax minMaxDistanceToQuadrant(const Geodesic &from, const Geodesic &toOrigin, const Geodesic &toEnd) {
        Geodesic fromRads{radsFromDegrees(from.latitude), radsFromDegrees(from.longitude)};
        Geodesic toOriginRads{radsFromDegrees(toOrigin.latitude), radsFromDegrees(toOrigin.longitude)};
        double toEndLat = std::min(toEnd.latitude, kMaxDegreesLatitude);
        double toEndLng = std::min(toEnd.longitude, kMaxDegreesLongitude);
        Geodesic toEndRads{radsFromDegrees(toEndLat), radsFromDegrees(toEndLng)};
        auto res = minMaxDistanceToQuadrant(fromRads.latitude, fromRads.longitude,
                                            toOriginRads.latitude, toEndRads.latitude,
                                            toOriginRads.longitude, toEndRads.longitude);
        return res;
    }

    static double distanceOnEarth(const Geodesic &from, const Geodesic &to) {
        return kEarthRadius * haversineDistance(from, to);
    }

    static double haversineDistance(const Geodesic &from, const Geodesic &to) {
        return radsHaversineDistance(radsFromDegrees(from), radsFromDegrees(to));
    }

    static double radsHaversineDistance(Geodesic from, Geodesic to) {
        double latSin = std::sin((from.latitude - to.latitude) * 0.5);
        double longSin = std::sin((from.longitude - to.longitude) * 0.5);
        double asinArg = std::sqrt(latSin * latSin + std::cos(from.latitude) * std::cos(to.latitude) * longSin * longSin);
        return 2 * std::asin(clamp(asinArg));
    }

    static double squareDegreeArea(Geodesic geodesic) {
        Geodesic topLeft{geodesic.latitude + 1, geodesic.longitude};
        Geodesic bottomRight{geodesic.latitude, geodesic.longitude + 1};
        return distanceOnEarth(geodesic, topLeft) * distanceOnEarth(geodesic, bottomRight);
    }

    static Distance distanceToEdgeNew(const Geodesic &from, const std::array<Geodesic, 2> &to, bool normalize) {
        assert(to[0].latitude < to[1].latitude && to[0].longitude < to[1].longitude);
        std::array<Distance, 4> distances;
        for (int i = 0; i < 2; ++i) { // i -- const coordinate
            for (int j = 0; j < 2; ++j) { // j -- value of i
                Geodesic sideFrom;
                Geodesic sideTo;
                sideFrom.c[i] = to[j].c[i];
                sideTo.c[i] = to[j].c[i];
                sideFrom.c[1-i] = to[0].c[1-i];
                sideTo.c[1-i] = to[1].c[1-i];
                assert(sideFrom.latitude <= sideTo.latitude && sideFrom.longitude <= sideTo.longitude);
                if (normalize) {
                    sideFrom.normalize();
                    sideTo.normalize();
                }
                assert(sideFrom.latitude <= sideTo.latitude && sideFrom.longitude <= sideTo.longitude);
                distances[i * 2 + j] = minMaxDistanceToQuadrant(from, sideFrom, sideTo).min;
            }
        }
        return *std::min_element(distances.begin(), distances.end());
    }

    static Distance distanceToEdge(const Geodesic &from, const Geodesic &toOrigin, const Geodesic &toEnd) {
        Distance toLeft = minMaxDistanceToQuadrant(from, {toOrigin.latitude, toOrigin.longitude}, {toEnd.latitude, toOrigin.longitude}).min;
        Distance toBottom = minMaxDistanceToQuadrant(from, {toOrigin.latitude, toOrigin.longitude}, {toOrigin.latitude, toEnd.longitude}).min;
        Distance toRight = minMaxDistanceToQuadrant(from, {toOrigin.latitude, toEnd.longitude}, {toEnd.latitude, toEnd.longitude}).min;
        Distance toTop = minMaxDistanceToQuadrant(from, {toEnd.latitude, toOrigin.longitude}, {toEnd.latitude, toEnd.longitude}).min;
        return std::min({toLeft, toBottom, toRight, toTop});
    }
};