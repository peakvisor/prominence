#ifndef TESTS_PROMINENCE_CALCULATION_H
#define TESTS_PROMINENCE_CALCULATION_H

#include "common.h"
#include "index.h"
#include "field_view.h"
#include "logging.hpp"
#include "dem_filename_formatter.h"
#include "dem_grid.h"
#include "qtree.h"
#include "frontier_queue.h"
#include "peak.h"
#include "peaks_collector.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <filesystem>

template <typename View>
void discreteRefine(const DEMGrid<View> &region, std::deque<Peak> &peaks) {
    Distance maxShift = 0;
    for (auto &peak : peaks) {
        if (peak.prominence <= 0) {
            continue;
        }

        auto planar = region.geodesicToPlanar(peak.geo());
        SpatialIndex spatial{planar, region.get(planar)};
        SpatialIndex refinedSpatial = spatial;
        SpatialIndex nextSpatial = region.getHighestNeighbour(spatial);
        while (refinedSpatial < nextSpatial) {
            std::swap(refinedSpatial, nextSpatial);
            nextSpatial = region.getHighestNeighbour(refinedSpatial);
        }

        if (spatial < refinedSpatial) {
            Geodesic newGeo = region.planarToGeodesic(refinedSpatial.planar);
            auto shift = DistanceTools::distanceOnEarth(newGeo, peak.geo());
            if (shift > maxShift) {
                maxShift = shift;
            }
            std::cout << "shifting " << shift << endl << peak.toDescriptiveString() << endl;
            peak.latitude = newGeo.latitude;
            peak.longitude = newGeo.longitude;
            if (peak.altitude != refinedSpatial.altitude) {
                std::cout << "rising\n";
            }
            peak.altitude = refinedSpatial.altitude;
            std::cout << "to\n" << peak.toDescriptiveString() << endl;
        }
    }
    std::cout << OUTSTR(maxShift) << endl;
}

template <typename View>
std::deque<Peak> loadPeaksForForRegion(const DEMGrid<View> &region, const std::string &inputDatabase) {
    return loadPeaks(inputDatabase, [&region](const Peak& peak) {
        return region.geodesicIsInside(peak.latitude, peak.longitude);
    });
}

void discreteRefine(const std::string &demsDir, const std::string &inputDatabase, const std::string &outputDatabase, const PlanarIndex from, const PlanarIndex to) {
    TiledDEMGrid region(demsDir, from, to, {kBigTiledFileSizeX, kBigTiledFileSizeY}, kHDDemSize);
    auto peaks = loadPeaksForForRegion(region, inputDatabase);
    discreteRefine(region, peaks);
    savePeaks(peaks, outputDatabase);
}

void printInfoIfNeeded(const uint64_t &visited, const uint64_t &totalPoints,
                       const Geodesic &geo, const double &visitedFromAbove,
                       const size_t &qtreeNodes, const Altitude &alt,
                       const Altitude &waterLevel) {

    static const size_t secondsPerInfo = 10;

    static auto startTime = clock();

    auto nextTime = clock();
    static auto previousTime = startTime;
    if (secondsBetween(previousTime, nextTime) < secondsPerInfo) {
        return;
    }
    previousTime = nextTime;
    auto timeSpent = secondsBetween(startTime, nextTime);
    auto visitedShareEstimate = double(visited) / totalPoints;
    flog() << std::setprecision(4) << 100 * visitedShareEstimate << "% completed in "
           << timeSpent << "s est time left: " << (timeSpent / visitedShareEstimate - timeSpent) / 60.
           << "m visiting " << geo
           << " a: " << alt
           << " fromAbove: " << (100. * visitedFromAbove) / visited
           << "%" << " qtree nodes: " << qtreeNodes
           << " waterLevel: " << waterLevel << std::endl;
}

bool outputBuggedVisit(const uint64_t &visited, const uint64_t &totalPoints,
                       const Seed &nextPoint) {
    if (visited == totalPoints + 1) {
        flog() << "Setting output zero for PlanarIndex to: " << nextPoint.planar() << std::endl;
        PlanarIndex::outputZero = nextPoint.planar();
    }
    flog() << "visited: " << visited << " visiting: " << nextPoint.planar() << " a: " << nextPoint.altitude() << std::endl;
    static const uint64_t kPointsToOutput = 50000;
    return visited < totalPoints + kPointsToOutput;
}

template <typename View>
std::deque<Peak> calculateMultiseedProminences(const DEMGrid<View> &region,
                                                   std::deque<Peak> &&inputPeaks,
                                                   const std::string &outputPeaksFile,
                                                   const Altitude &oceanFloor = kOceanFloor,
                                                   const Altitude &additionalPeaksCutoff = 100) {
    bool mainProcedure = oceanFloor == kOceanFloor;
    bool searchingForAdditionalPeaks = additionalPeaksCutoff > 0;
    auto startTime = clock();
    PeaksContainer peaksContainer(region.converter, std::move(inputPeaks), outputPeaksFile);
    peaksContainer.init(region);

    FrontierQueue frontier(peaksContainer.seeds, oceanFloor);
    auto frontierPusher = [&region, &frontier](const PlanarIndex planar) -> bool {
        Altitude altitude = region.get(planar);
        return frontier.push(planar, altitude);
    };

    QTree qtree(region.indexFrom(), region.indexTo(), peaksContainer, region.converter, frontierPusher);

    PeaksCollector peaksCollector(additionalPeaksCutoff);

    SpatialIndex lastSpatialNotFromAbove{{0, 0}, kNoAltitude};
    const uint64_t totalPoints = uint64_t(region.filesCount) * region.cellsResolution.x * region.cellsResolution.y + peaksContainer.seeds.size();

    uint64_t visited = 0;
    uint64_t visitedFromAbove = 0;

    auto waterLevel = frontier.getWaterLevel();
    auto previousWaterLevel = waterLevel;

    while ((waterLevel = frontier.getWaterLevel()) > oceanFloor) {
        assert(previousWaterLevel >= waterLevel);
        previousWaterLevel = waterLevel;
        const auto &nextPoint = frontier.pop();

        if (nextPoint.altitude() <= oceanFloor) {
            flog() << "finishing due to next point being lower that ocean floor\n";
            break;
        }

        if (peaksContainer.finished()) {
            flog() << "finishing because peaksContainer is finished\n";
            break;
        }

        if (nextPoint.islandOwner != kNoIsland && nextPoint.islandOwner != 0 && mainProcedure) {
            peaksContainer.setIsolation(nextPoint.islandOwner, qtree.computeIsolation(nextPoint.planar()));
        }

        assert(nextPoint.altitude() >= waterLevel);

        ++visited;
        static const size_t kVisitsPerCheckTime = 10000;
        if (visited % kVisitsPerCheckTime == 0) {
            printInfoIfNeeded(visited, totalPoints, region.planarToGeodesic(nextPoint.planar()),
                visitedFromAbove, qtree.nodesSize(), nextPoint.altitude(), waterLevel);
        }

        if (visited > totalPoints) {
            bool shouldContinueGatheringInfo = outputBuggedVisit(visited, totalPoints, nextPoint);
            if (!shouldContinueGatheringInfo) {
                peaksContainer.markAllLeftUnreliable();
                flog() << "err: finishing due to visit number exceding totalPoints\n";
                break;
            }
        }

        auto island = qtree.visitPoint(nextPoint);
        auto owner = peaksContainer.getIslandOwner(island);
        assert(island == owner);

        if (nextPoint.altitude() > waterLevel && searchingForAdditionalPeaks) {
            ++visitedFromAbove;
            if (peaksContainer.hasReliableProminence(island)) {
                if (peaksCollector.inited == false) {
                    peaksCollector.init(lastSpatialNotFromAbove, nextPoint.spatial, island);
                    assert(lastSpatialNotFromAbove.altitude == waterLevel);
                    if (peaksContainer.getPeak(island).altitude < waterLevel) {
                        std::cout << peaksContainer.getPeak(island).toDescriptiveString() << std::endl;
                        std::cout << waterLevel << std::endl;
                        assert(false);
                    }
                }
                assert(peaksCollector.parent == island);
                peaksCollector.visit(region, nextPoint.spatial);
            } else {
                assert(!peaksCollector.inited || peaksCollector.parent == island);
                peaksCollector.peaks.clear();
            }
        } else if (searchingForAdditionalPeaks) {
            if (peaksCollector.inited) {
                auto &parentPeak = peaksContainer.getPeak(peaksCollector.parent);
                if (peaksCollector.peaks.empty() && SpatialIndex{region.geodesicToPlanar(parentPeak.geo()), parentPeak.altitude} < peaksCollector.highest) {
                    peaksCollector.addPeak(region, peaksCollector.highest);
                }

                if (!peaksCollector.peaks.empty()) {
                    std::deque<Peak> filteredPeaks;

                    if (peaksCollector.peaks.size() > 1) {
                        FileLog::setPrefix("extra ");
                        auto notFilteredPeaks = calculateMultiseedProminences(region, std::move(peaksCollector.peaks), "tmp_extra_peaks.txt", peaksCollector.keycol.altitude);

                        size_t highestIndex = kNoIndexSizeT;
                        for (size_t i = 0; i < notFilteredPeaks.size(); ++i) {
                            if (highestIndex == kNoIndexSizeT || notFilteredPeaks[highestIndex].definitiveAltitudeLess(notFilteredPeaks[i])) {
                                highestIndex = i;
                            }
                        }

                        for (const auto &peak : notFilteredPeaks) {
                            if (peak.prominence >= additionalPeaksCutoff || peak.id == notFilteredPeaks.at(highestIndex).id) {
                                filteredPeaks.push_back(peak);
                            }
                        }
                        notFilteredPeaks.clear();
                        FileLog::setPrefix("");
                    } else {
                        filteredPeaks.push_back(peaksCollector.peaks.front());
                        filteredPeaks.front().islandParent = filteredPeaks.front().id;
                    }

                    for (const auto &peak : filteredPeaks) {
                        assert(peak.altitude > peaksCollector.keycol.altitude);
                    }

                    peaksContainer.addExtraPeaks(filteredPeaks, peaksCollector.keycol, peaksCollector.parent);

                }
                peaksCollector.deinit();
            }

            lastSpatialNotFromAbove = nextPoint.spatial;
        }

        frontier.liftWater();

        if (peaksContainer.failed) {
            break;
        }
    }

    peaksContainer.postProcess();

    if (mainProcedure) {
        flog() << "total visited: " << visited << " fromAbove: " << visitedFromAbove << std::endl;
        peaksContainer.d__outputProminenceMapping();
        logTime("Total", startTime);
        flog() << "output peaks at " << outputPeaksFile << std::endl;
    }

    return peaksContainer.unloadPeaks();
}

template <typename View>
void calculateProminences(const std::string &demsDirectory,
                          const std::string &inputPeaksFile,
                          const std::string &outputPeaksFile,
                          const PlanarIndex &from, const PlanarIndex &to,
                          PlanarIndex fileSize, Index cellsResolution,
                          const Altitude &extraPeaksProminenceCutoff = 100) {
    DEMGrid<View> region{demsDirectory, from, to, fileSize, cellsResolution, true};
    auto peaks = loadPeaks(inputPeaksFile, [&region](const Peak& peak) {
        return region.geodesicIsInside(peak.latitude, peak.longitude);
    });
    assert(peaks.size() > 0);
    std::cout << "total seeds: " << peaks.size() << std::endl;
//    getIdToIndexMap(peaks);
    auto resultPeaks = calculateMultiseedProminences<View>(region, std::move(peaks), outputPeaksFile);
}

template <typename DEM>
std::deque<Peak> calculateProminences(const DEM &dem, const std::filesystem::path &inputPeaksFile, const std::string &backupOutputPeaksFile) {
    auto peaks = loadPeaks(inputPeaksFile, [&dem](const Peak& peak) {
        return dem.geodesicIsInside(peak.latitude, peak.longitude);
    });
    return calculateMultiseedProminences(dem, std::move(peaks), backupOutputPeaksFile);
}

template <typename DEM>
std::deque<Peak> calculateProminences(const DEM &dem, std::deque<Peak> &peaks, const std::string &backupOutputPeaksFile) {
    return calculateMultiseedProminences(dem, std::move(peaks), backupOutputPeaksFile);
}

int applyProminenceCalculation(int argc, char *argv[]) {
    LogGuard log("prominence calculation");

    std::string mode = argc > 1 ? std::string(argv[1]) : "in_code"s;
    if (mode == "pc") {
        bool usingHDDems = *argv[2] == 'h';
        std::string demDirectory = std::string(argv[3]) + "/";
        PlanarIndex from = {std::atoi(argv[4]), std::atoi(argv[5])};
        PlanarIndex to = argc > 7 ? PlanarIndex{std::atoi(argv[6]), std::atoi(argv[7])} : from;
        std::string inputPeaks = argc > 8 ? std::string(argv[8]) : "peaks.txt";
        std::string outputPeaks = argc > 9 ? std::string(argv[9]) : "output_peaks.txt";
        std::string logFile = argc > 10 ? std::string(argv[10]) : "pcal.log";

        FileLog::setUpFileLog(logFile);

        PlanarIndex fileSize = usingHDDems ? PlanarIndex{kBigTiledFileSizeX, kBigTiledFileSizeY} : PlanarIndex{1200, 1200};
        Index cellsResolution = usingHDDems ? 3600 : 1200;

        if (usingHDDems) {
            calculateProminences<TiledFieldView>(demDirectory, inputPeaks, outputPeaks, from, to, fileSize, cellsResolution);
        } else {
            calculateProminences<FieldView>(demDirectory, inputPeaks, outputPeaks, from, to, fileSize, cellsResolution);
        }
    } else if (mode == "drefine") {
        discreteRefine(argv[2], argv[3], argv[4], {-47, 166}, {-35, 178});
    } else if (mode == "in_code") {
        std::string demDirectory = "/Users/ka/pc/prod/hd_dems/";
        std::string inputPeaks = "/Users/ka/pc/prod/pipeline/2020-12-22/refined_peaks.txt";
        std::string outputPeaks = "/Users/ka/pc/prod/tmp/output_peaks.txt";
        FileLog::setUpFileLog("/Users/ka/pc/prod/tmp/log.txt");
        PlanarIndex from{33, 126};
        PlanarIndex to{33, 126};
        calculateProminences<TiledFieldView>(demDirectory, inputPeaks, outputPeaks, from, to, {kBigTiledFileSizeX, kBigTiledFileSizeY}, 3600);
    }
    return 0;
}

#endif //TESTS_PROMINENCE_CALCULATION_H
