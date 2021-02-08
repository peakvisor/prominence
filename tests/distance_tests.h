#ifndef TESTS_DISTANCE_TESTS_H
#define TESTS_DISTANCE_TESTS_H
#include <random>
#include <iostream>

constexpr double operator"" _deg(long double d) { return radsFromDegrees(d); }

double randomLat(double min = -M_PI_2, double max = M_PI_2) {
    static std::mt19937 generator(std::random_device{}());
    return std::uniform_real_distribution<double>(min, max)(generator);
}

double randomLng(double min = -M_PI, double max = M_PI) {
    static std::mt19937 generator(std::random_device{}());
    return std::uniform_real_distribution<double>(min, max)(generator);
}

void testDistanceToSegmentOnMeridianOneTime(double fromLat, double fromLng, double meridian, double toBottomLat, double toTopLat, bool isMinDistance) {
    bool failed = false;
    try {
        auto dist = DistanceTools::distanceToSegmentOnMeridian(fromLat, fromLng, meridian, toBottomLat, toTopLat,
            isMinDistance);
        Geodesic segmentCenter = {(toBottomLat + toTopLat) / 2, meridian};
        auto distToCenter = DistanceTools::radsHaversineDistance({fromLat, fromLng}, segmentCenter);
        bool failedComparison = isMinDistance ? dist > distToCenter + DistanceTools::kRadEps : dist < dist - DistanceTools::kRadEps;
        if (failedComparison) {
            failed = true;
            std::cout << OUTSTR(isMinDistance) << std::endl;
            std::cout << OUTSTR(dist) << " " << OUTSTR(distToCenter) << std::endl;
        }
    } catch (const std::logic_error &e) {
        std::cout << e.what() << std::endl;
        failed = true;
    }

    if (failed) {
        std::cout << fromLat << ", " << fromLng << ", " << meridian << ", "
                  << toBottomLat << ", " << toTopLat << ", " << isMinDistance << std::endl;
        std::cout << degreesFromRads(fromLat) << " " << degreesFromRads(fromLng) << " meridian: " << degreesFromRads(meridian) << " segment: ["
                  << degreesFromRads(toBottomLat) << ", " << degreesFromRads(toTopLat) << "]\n";
        throw std::logic_error("distanceToSegmentOnMeridian failed");
    }
}

void testDistance() {
    int iterations = 100000;
    double step = 0.01;
    double currentPercent = step;
    for (int i = 0; i <  iterations; ++i) {
        if (i > currentPercent * iterations) {
//            std::cout << currentPercent * 100 << "%\n";
            currentPercent += step;
        }
        auto fromLat = randomLat();
        auto fromLng = randomLng();
        auto meridian = randomLng();
        auto toBottomLat = randomLat();
        auto toTopLat = randomLat();

        if (toBottomLat > toTopLat) {
            std::swap(toTopLat, toBottomLat);
        }
        auto isMinDistance = bool(i % 2);
        try {
            testDistanceToSegmentOnMeridianOneTime(fromLat, fromLng, meridian, toBottomLat, toTopLat, isMinDistance);
        } catch (const std::logic_error &e) {
            std::cout << "failed iteration: " << i << std::endl;
            break;
        }
    }
}

#endif //TESTS_DISTANCE_TESTS_H
