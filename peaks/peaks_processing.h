#ifndef TESTS_PEAK_LIST_TOOLS_H
#define TESTS_PEAK_LIST_TOOLS_H

#include "index.h"
#include "distance.h"
#include "peak.h"
#include "logging.hpp"
#include "peak_tree.hpp"
#include "passes_processing.h"
#include "converter.h"

#include <random>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <map>
#include <algorithm>
#include <set>
#include <filesystem>

void applyRefine(const std::string &inputDatabase, const std::string &refinementsDatabase,
        const std::string &outputDatanase) // TODO: test
{
    auto peaks = loadPeaks(inputDatabase);
    auto refs = loadPeaks(refinementsDatabase);

    updateSortedPeaks(peaks, refs, [](Peak &peak, const Peak &ref) {
        peak.latitude = ref.latitude;
        peak.longitude = ref.longitude;
        peak.altitude = ref.altitude;
    });
    savePeaks(peaks, outputDatanase);
}

void checkPostProcessed(const Path &pipelineDir) {
    std::string postProcessedDatabase = pipelineDir / "post_processed.txt";
    std::string rawDatabase = pipelineDir / "input_peaks.txt";
    std::string vipIdsList = pipelineDir / "vip_ids.txt";
    std::string fixedAltitudeList = pipelineDir / "fixed_alt.txt";
    auto peaks = loadPeaks(postProcessedDatabase);
    auto vipIds = loadPeaks(vipIdsList);
    auto lookForOnlySecond = [&](size_t fi, size_t si, DualState state) {
        if (state == kDualStateSecond) {
            std::cout << "vip id missed: " << vipIds.at(si).id << std::endl;
        }
    };
    dualIteration(peaks, vipIds, lookForOnlySecond);
    auto refPeaks = loadPeaks(rawDatabase);
    auto fixedAltList = loadPeaks(fixedAltitudeList);
    std::deque<Peak> fixedAltPeaks;

    auto fillSavedPeaks = [&](size_t fi, size_t si, DualState state) {
        if (state == kDualStateBoth) {
            fixedAltPeaks.emplace_back(refPeaks.at(fi));
        }
    };

    dualIteration(refPeaks, fixedAltList, fillSavedPeaks);

    auto checkFixedAlt = [&](size_t fi, size_t si, DualState state) {
        if (state == kDualStateBoth) {
            auto &peak = peaks.at(fi);
            auto &refPeak = fixedAltPeaks.at(si);
            Altitude diff = peak.altitude - refPeak.altitude;
            if (diff != 0) {
                if (fixedAltPeaks.at(si).altitude <= 0) {
                    std::cout << refPeak.toDescriptiveString() << " has zero altitude in input\n";
                } else {
                    std::cout << peak.toDescriptiveString();
                    std::cout << "\tshould be " << refPeak.altitude << " high, diff: " << diff
                              << std::endl;
                }
            }
        }
    };

    dualIteration(peaks, fixedAltPeaks, checkFixedAlt);

    PeakTree tree(peaks);
    for (auto i = 0u; i < peaks.size(); ++i) {
        auto peak = peaks.at(i);
        if (!peak.isExtra() && peak.islandParent != kNoId) {
            if (peaks.at(tree.nodes.at(i).parent).isExtra()) {
                std::cout << peak.toDescriptiveString() << " has extra parent\n";
            }
        }
    }

    for (const auto &peak : peaks) {
        if (peak.prominence > 0) {
            if (peak.ilpLatitude == kNoDegrees) {
                std::cout << peak.toDescriptiveString() << " without isolation\n";
            }
        }
    }
}

void printTopByCallback(std::deque<Peak> &peaks, std::function<double(const Peak&)> callback, size_t top) {
    std::sort(peaks.begin(), peaks.end(), [&callback](const Peak &f, const Peak &s) {
        return callback(f) > callback(s);
    });
    for (size_t i = 0; i < top; ++i) {
        std::cout << i + 1 << "\t" << peaks[i].toDescriptiveString() << std::endl;
    }
}

void peaksNear(const std::string &inputDatabse, Geodesic co, uint8_t peaksCount = 10) {
    auto peaks = loadPeaks(inputDatabse);
    std::deque<std::pair<Distance, size_t>> nearestPeaks;
    std::cout << "looking for peaks near " << co << endl;
    for (size_t peakIndex = 0; peakIndex < peaks.size(); ++peakIndex) {
        const auto &peak = peaks[peakIndex];
        auto dist = DistanceTools::distanceOnEarth(co, peak.geo());
        for (size_t nearestPeakIndex = 1u; nearestPeakIndex < nearestPeaks.size(); ++nearestPeakIndex) {
            assert(nearestPeaks[nearestPeakIndex - 1].first < nearestPeaks[nearestPeakIndex].first);
        }
        if (nearestPeaks.size() <= peaksCount || nearestPeaks.back().first > dist) {
            nearestPeaks.emplace_back(dist, peakIndex);
            std::sort(nearestPeaks.begin(), nearestPeaks.end());
            if (nearestPeaks.size() > peaksCount) {
                nearestPeaks.pop_back();
            }
        }
    }
    for (const auto &distPeak : nearestPeaks) {
        std::cout << distPeak.first << ":\t" << peaks.at(distPeak.second).toDescriptiveString() << endl;
    }
    std::cout << "finish peaks near\n";
}

void findProminentChanges(const std::string &oldDatabase, const std::string &newDatabase) {
    static Distance kDistanceCutoff = 0.5;
    static Altitude kAltitudeDiffCutoff = 10;
    static double kAltitudePropDiffCutoff = 0.2;
    static Altitude kProminenceDiffCutoff = 10;
    static double kProminencePropDiffCutoff = 0.2;
    std::cout << "diff margins:\n" << "kDistanceCutoff = " << kDistanceCutoff << " kAltitudeDiffCutoff = " << kAltitudeDiffCutoff
              << " kAltitudeDiffCutoff = " << kAltitudeDiffCutoff << " kProminenceDiffCutoff = " << kProminenceDiffCutoff
              << " kProminencePropDiffCutoff =  " << kProminencePropDiffCutoff << endl;
    std::cout << "diff is counted if loweredProminence && (areDistant || bigAltDiff || bigPromDiff)" << endl;
    std::cout << "bigAltDiff = altitudePropDiff > kAltitudePropDiffCutoff && altitudeDiff > kAltitudeDiffCutoff\n";
    std::cout << "bigPromDiff = prominencePropDiff > kProminencePropDiffCutoff && prominenceDiff > kProminenceDiffCutof\n";
    auto oldPeaks = loadPeaks(oldDatabase);
    auto newPeaks = loadPeaks(newDatabase);
    sortById(oldPeaks);
    sortById(newPeaks);
    std::vector<std::pair<size_t, size_t>> diffs;
    auto compareCallback = [&](size_t oldIndex, size_t newIndex, DualState dualState) {
        if (dualState == kDualStateFirst) {
            return;
        }
        if (newIndex >= newPeaks.size()) {
            std::cout << newIndex << endl;
            std::cout << newPeaks.size() << endl;
            assert(false);
        }
        auto &newPeak = newPeaks.at(newIndex);
        bool shouldPost = false;

        if (dualState == kDualStateBoth) {
            assert(oldIndex < oldPeaks.size());
            auto &oldPeak = oldPeaks.at(oldIndex);
            bool loweredProminence = newPeak.prominence < oldPeak.prominence && oldPeak.prominence != 0 && oldPeak.prominence != kNoAltitude;
            if (loweredProminence) {
                double distance = DistanceTools::distanceOnEarth(oldPeak.geo(), newPeak.geo());
                double altitudePropDiff = std::abs(static_cast<double>(oldPeak.altitude + 1.) / (newPeak.altitude + 1.) - 1);
                Altitude altitudeDiff = std::abs(newPeak.altitude - oldPeak.altitude);
                double prominencePropDiff = std::abs(static_cast<double>(oldPeak.prominence + 1.) / (newPeak.prominence + 1.) - 1);
                Altitude prominenceDiff = std::abs(newPeak.prominence - oldPeak.prominence);

                bool areDistant = distance > kDistanceCutoff;
                bool bigAltDiff = altitudePropDiff > kAltitudePropDiffCutoff && altitudeDiff > kAltitudeDiffCutoff;
                bool bigPromDiff = prominencePropDiff > kProminencePropDiffCutoff && prominenceDiff > kProminenceDiffCutoff;

                if (areDistant || bigAltDiff || bigPromDiff) {
                    shouldPost = true;
                }
            }

            if (shouldPost) {
                diffs.emplace_back(oldIndex, newIndex);
            }
        }
    };

    dualIteration(oldPeaks, newPeaks, compareCallback, true);
    std::sort(diffs.begin(), diffs.end(), [&oldPeaks](const auto &lhs, const auto &rhs) {
        return oldPeaks.at(lhs.first).prominence > oldPeaks.at(rhs.first).prominence;
    });

    size_t count = 0;
    for (const auto &diff : diffs) {
        auto oldIndex = diff.first;
        auto newIndex = diff.second;
        auto &oldPeak = oldPeaks.at(oldIndex);
        auto &newPeak = newPeaks.at(newIndex);
        std::cout << "diff " << count++ << endl;
        std::cout << "distance: " << DistanceTools::distanceOnEarth(oldPeak.geo(), newPeak.geo()) << " prom diff: " << oldPeak.prominence - newPeak.prominence << endl;
        std::cout << oldPeak.toDescriptiveString() << endl;
        std::cout << newPeak.toDescriptiveString() << endl;
    }
    std::cout << "total diffs: " << diffs.size() << endl;
    std::cout << "finish comparing datasets\n";

}

void compareDatabases(const std::string &firstDatabase, const std::string &secondDatabase) {
    auto firstPeaks = loadPeaks(firstDatabase);
    auto secondPeaks = loadPeaks(secondDatabase);

    dualIteration(firstPeaks, secondPeaks, [&firstPeaks, &secondPeaks](size_t fi, size_t si, DualState dualState) {
        switch (dualState) {
            case kDualStateFirst:
                std::cout << "1st only: " << firstPeaks.at(fi).toDescriptiveString() << endl;
                break;
            case kDualStateSecond:
                std::cout << "2nd only: " << secondPeaks.at(si).toDescriptiveString() << endl;
                break;
            case kDualStateBoth:
            default:
                break;
        }
    });
}





void standardisePeaks(std::deque<Peak> &peaks) {
    PeakTree tree(peaks);
    tree.resolveSameAltitudeParenting(false);
    std::unordered_map<PeakId, PeakId> unnamingMap;
    std::unordered_set<PeakId> newIds;
    Converter converter{{-90, -180}, 3600};
    for (const auto &peak : peaks) {
        if (peak.prominence == 0 || peak.altitude == 0) {
            continue;
        }
        PlanarIndex planar = converter.geodesicToPlanar(peak.geo());
        Geodesic geo = converter.planarToGeodesic(planar);
        PeakId newId = extraPeakId(geo);
        auto ins = unnamingMap.try_emplace(peak.id, newId);
        assert(ins.second);
        auto newIdsIns = newIds.emplace(newId);
        if (!newIdsIns.second) {
            std::cout << peak.toDescriptiveString() << std::endl;
            std::cout << *newIdsIns.first << std::endl;
            assert(newIdsIns.second);
        }
    }

    std::deque<Peak> standardizedPeaks;

    for (auto &peak : peaks) {
        if (peak.prominence <= 0) {
            continue;
        }
        auto it = unnamingMap.find(peak.id);
        if (it == unnamingMap.end()) {
            std::cout << "no " << OUTSTR(peak.id) << " in unnaming map\n";
            std::cout << peak.toDescriptiveString() << std::endl;
        } else {
            peak.id = it->second;
        }
        if (peak.islandParent != kNoId) {
            auto pit = unnamingMap.find(peak.islandParent);
            if (pit == unnamingMap.end()) {
                std::cout << "no " << OUTSTR(peak.islandParent) << " in unnaming map\n";
                std::cout << peak.toDescriptiveString() << std::endl;
            } else {
                peak.islandParent = pit->second;
            }
        }

        standardizedPeaks.push_back(peak);
    }
    std::exchange(peaks, standardizedPeaks);
}

void filterByCoordinates(const std::string &inputPeaksFile, const std::string &outputPeaksFile, const Geodesic &from, const Geodesic &to) {
    auto peaks = loadPeaks(inputPeaksFile, from, to);
    savePeaks(peaks, outputPeaksFile);
}

void applyPeakListTools(int argc, char *argv[]) {
    auto startTime = clock();

    if (argc > 1) {
        std::string mode = argv[1];
        if (mode == "apply_refine") {
            std::string inputDatabase = std::string(argv[2]);
            std::string refsDatabase = std::string(argv[3]);
            std::string outputDatabase = std::string(argv[4]);
            applyRefine(inputDatabase, refsDatabase, outputDatabase);
        } else if (mode == "peaks_near") {
            if (argc == 6) {
                peaksNear(argv[2], {std::stod(argv[3]), std::stod(argv[4])}, std::atoi(argv[5]));
            } else {
                std::cout << "peaksNear(const std::string &inputDatabse, Geodesic co, uint8_t peaksCount)\n";
            }
        } else if (mode == "find_prominent_changes") {
            std::cout << "find_prominent_changes\n";
            findProminentChanges(argv[2], argv[3]);
        } else if (mode == "compare_databases") {
            std::cout << "compare databases\n";
            compareDatabases(argv[2], argv[3]);
        } else if (mode == "top_prominent") {
            auto peaks = loadPeaks(argv[2]);
            size_t count = std::stoi(argv[3]);
            printTopByCallback(peaks, [](const Peak &peak) { return peak.prominence; }, count);
        } else if (mode == "top_prominent_extras") {
            auto peaks = loadPeaks(argv[2]);
            size_t count = std::stoi(argv[3]);
            printTopByCallback(peaks, [](const Peak &peak) { return peak.prominence * peak.isExtra(); }, count);
        } else if (mode == "standardize") {
            auto peaks = loadPeaks(argv[2]);
            standardisePeaks(peaks);
            savePeaks(peaks, argv[3]);
        } else if (mode == "check_post_processed") {
            checkPostProcessed(argv[2]);
        } else if (mode == "distance") {
            std::cout << DistanceTools::distanceOnEarth({std::stod(argv[2]), std::stod(argv[3])},
                    {std::stod(argv[4]), std::stod(argv[5])}) << std::endl;
        } else if (mode == "filter_by_coordinates") {
            filterByCoordinates(argv[2], argv[3], {std::stod(argv[4]), std::stod(argv[5])}, {std::stod(argv[6]), std::stod(argv[7])});
        }
        else if (mode == "in_code") {
            std::cout << "in_code mode\n";
        } else {
            std::cout << "wrong mode: " << mode << endl;
        }
    } else {
        std::cout << "missing mode\n";
    }

    std::cout << secondsBetween(startTime, clock()) << "s to finish\n";
}

#endif //TESTS_PEAK_LIST_TOOLS_H
