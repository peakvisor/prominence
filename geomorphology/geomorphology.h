//
// Created by Egor on 05.01.2021.
//

#ifndef TESTS_GEOMORPHOLOGY_H
#define TESTS_GEOMORPHOLOGY_H

#include "common.h"
#include "dem_grid.h"
#include "field_view.h"
#include "index.h"
#include "logging.hpp"
#include "peak.h"
#include "ridge_calculator.h"

#include <iostream>

void calcRidges(PlanarIndex from, PlanarIndex to, const std::string &demsDir, const std::string &referenceDatabase,
                const std::string &ridgesOutputDir, Index margin = 1) {
    logTime();
    auto regionFrom = from - kDiagonalShift * margin;
    auto regionTo = to + kDiagonalShift * margin;
    TiledDEMGrid region(demsDir, regionFrom, regionTo, {kBigTiledFileSizeX, kBigTiledFileSizeY}, 3600, true, false);
    logTime("region");

    RidgeCalculator<FieldView_<LinearIndexer, uint16_t>> calc(region, ridgesOutputDir, margin);

    logTime("setup calculator");

    calc.altitudeGetter = std::make_unique<AltitudeGetter>([&region](const PlanarIndex &planar) {
        return region.get(planar);
    });

    auto peaks = loadPeaks(referenceDatabase, [&region](const Peak& peak) {
        return region.geodesicIsInside(peak.latitude, peak.longitude);
    });

    logTime("peaks load");

    calc.addPeaksTree(peaks, region.converter);

    logTime("add peaks");

    calc.calculate();

    logTime("calculation");

    auto size = calc.peakNodes.size();

    std::cout << "total peaks nodes: " << size << endl;

}

void applyGeomorphologyTools(int argc, char * argv[]) {
    auto startTime = clock();
    std::string mode = argv[1];
    if (mode == "simple_ridge_calculator") {
        if (argc == 9) {
            calcRidges({std::stoi(argv[2]), std::stoi(argv[3])}, {std::stoi(argv[4]), std::stoi(argv[5])},
                argv[6], argv[7], argv[8], 1);
        } else {
            std::cout << "calcRidges(PlanarIndex from, PlanarIndex to, string demsDir, std::string inputDatabase, string ridgesOutputDir)\n";
        }
    } else if (mode == "in_code") {
        calcRidges({-47, 166}, {-35, 178},
            "/Users/ka/pc/prod/hd_dems/",
            "/Users/ka/pc/prod/peaks/nz_resolved_output_peaks.txt",
            "/Users/ka/pc/prod/nz_ridges/");
    } else {
        std::cout << "wrong mode " << mode << endl;
        std::cout << "modes: simple_ridge_calculator, in_code\n";
    }


    std::cout << "total time: " << secondsBetween(startTime, clock()) << "s\n";
}

#endif //TESTS_GEOMORPHOLOGY_H
