#ifndef PROMINENCE_CALCULATOR_CONVERTER_H
#define PROMINENCE_CALCULATOR_CONVERTER_H

#include "index.h"

struct Converter {
    Converter(PlanarIndex degreesShift_, PlanarIndex planarIndexResolution_)
        : degreesShift{degreesShift_},
        planarIndexResolution{planarIndexResolution_},
        latitudeWorldLength{180 * planarIndexResolution_.x},
        longitudeWorldLength{360 * planarIndexResolution_.y},
        minLatitude{(-90 + degreesShift_.latitude) * planarIndexResolution_.x},
        maxLatitude{minLatitude + latitudeWorldLength},
        minLongitude{(-180 + degreesShift_.longitude) * planarIndexResolution_.y},
        maxLongitude{minLongitude + longitudeWorldLength},
        inversePlanarIndexResolution{1. / planarIndexResolution_.x, 1. / planarIndexResolution_.y} {}

    Converter(PlanarIndex degreesShift_, Index planarIndexResolution_) : Converter(degreesShift_, PlanarIndex{planarIndexResolution_, planarIndexResolution_}) {}

    [[nodiscard]] inline Geodesic planarToGeodesic(PlanarIndex planar) const {
        return {planar.latitude * inversePlanarIndexResolution.x - degreesShift.latitude,
            planar.longitude * inversePlanarIndexResolution.x - degreesShift.longitude};
    }

    [[nodiscard]] inline PlanarIndex geodesicToPlanar(Geodesic geo) const {
        double shiftedLat = geo.latitude + degreesShift.latitude;
        double shiftedLong = geo.longitude + degreesShift.longitude;
        return {
            static_cast<int>(std::round(shiftedLat * planarIndexResolution.x)),
            static_cast<int>(std::round(shiftedLong * planarIndexResolution.y)),
        };
    }

    [[nodiscard]] inline double indexToDegreeArea(uint64_t indexArea) const {
        return indexArea * inversePlanarIndexResolution.x * inversePlanarIndexResolution.y;
    }

    [[nodiscard]] inline bool isNormalized(const PlanarIndex &planar) const {
        if (planar.latitude > maxLatitude) {
            return false;
        }

        if (planar.latitude < minLatitude) {
            return false;
        }

        if (planar.longitude >= maxLongitude) {
            return false;
        }

        if (planar.longitude < minLongitude) {
            return false;
        }
        return true;
    }

    inline void normalize(PlanarIndex &planar) const {
        if (planar.latitude > maxLatitude) {
            assert(planar.latitude < 2 * maxLatitude);
            planar.latitude = 2 * maxLatitude - planar.latitude;
            planar.longitude += 180;
        }

        if (planar.latitude < minLatitude) {
            planar.latitude = 2 * minLatitude - planar.longitude;
            planar.longitude += 180;
        }

        while (planar.longitude >= maxLongitude) {
            planar.longitude -= longitudeWorldLength;
        }

        while (planar.longitude < minLongitude) {
            planar.longitude += longitudeWorldLength;
        }

        assert(isNormalized(planar));
    }

    PlanarIndex convertToSignedUniversal(const PlanarIndex &planar) {
        return planar - degreesShift;
    }

    PlanarIndex convertFromSignedUniversal(const PlanarIndex &planar) {
        return planar + degreesShift;
    }

    const PlanarIndex degreesShift;
    const PlanarIndex planarIndexResolution;
    const Index latitudeWorldLength;
    const Index longitudeWorldLength;
    const Index minLatitude;
    const Index maxLatitude;
    const Index minLongitude;
    const Index maxLongitude;
    const Geodesic inversePlanarIndexResolution;
};

#endif //PROMINENCE_CALCULATOR_CONVERTER_H
