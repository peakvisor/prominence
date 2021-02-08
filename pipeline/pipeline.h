#ifndef PIPELINE_H
#define PIPELINE_H

#include "common.h"
#include "logging.hpp"
#include "peak_tree.hpp"
#include "peaks_processing.h"
#include "passes_processing.h"

#include <string>
#include <filesystem>
#include <functional>
#include <optional>

void removeZeroMarkedPeaks(std::deque<Peak> &peaks) {
    for (auto i = 0u; i < peaks.size(); ++i) {
        if (!peaks[i].id) {
            std::exchange(peaks[i], peaks.back());
            peaks.pop_back();
            --i;
        }
    }
}

struct ResourcesPaths {
    const Path resourcesDirPath;
    Path mounts() { return resourcesDirPath / "mounts.txt"; }
    Path vipIds() { return resourcesDirPath / "vip_ids.txt"; }
    Path fixedAltitude() { return resourcesDirPath / "fixed_alt.txt"; }
    Path passes() { return resourcesDirPath / "passes.txt"; }
};

struct Processor {
    virtual void process(std::deque<Peak> &peaks) = 0;
    virtual std::string getName() = 0;
    virtual ~Processor() = default;
};

struct RemoveDuplicateIdsProcessor : Processor {
    std::string getName() final { return "remove_duplicate_ids"; }

    void process(std::deque<Peak> &peaks) final {
        sortById(peaks);
        for (auto i = 0u; i < peaks.size() - 1; ++i) {
            auto &peak = peaks.at(i);
            auto &nextPeak = peaks.at(i + 1);
            if (peak.id == nextPeak.id) {
                if (peak.isExtra() == nextPeak.isExtra()) {
                    throw std::logic_error("same id and same extraseness peaks");
                }
                if (peak.isExtra()) {
                    peak.name = nextPeak.name;
                }
                if (nextPeak.isExtra()) {
                    nextPeak.name = peak.name;
                }
                if (int(peak.prominence == 0) + int(nextPeak.prominence == 0) != 1) {
                    throw std::invalid_argument("nontribial doubles in removeDuplicateIds");
                }
                auto indexToRemove = peak.prominence == 0 ? i : i + 1;
                peaks.erase(peaks.begin() + indexToRemove);
            }
        }
    }
};

struct AddMissingPeaksProcessor : Processor {
    Path referencePeaksFile;
    explicit AddMissingPeaksProcessor(Path referencePeaksFile) : referencePeaksFile(std::move(referencePeaksFile)) {}

    std::string getName() final { return "add_missing"; }

    void process(std::deque<Peak> &peaks) final {
        auto referencePeaks = loadPeaks(referencePeaksFile);
        std::deque<Peak> missedPeaks;
        auto lambda = [&missedPeaks, &referencePeaks](size_t fi, size_t si, DualState state) {
            if (state == kDualStateSecond) {
                missedPeaks.push_back(referencePeaks[si]);
            }
        };
        dualIteration(peaks, referencePeaks, lambda);
        peaks.insert(peaks.end(), missedPeaks.begin(), missedPeaks.end());
    }
};

std::deque<std::string> readStringsFromFile(const std::string &stringFile) {
    std::deque<std::string> res;
    std::ifstream input{stringFile};
    std::string nextLine;
    while (std::getline(input, nextLine)) {
        res.push_back(nextLine);
    }

    return res;
}

std::deque<std::string> split(const std::string &s) {
    std::deque<std::string> res;
    std::istringstream scanner;
    std::string nextComponent;
    scanner.str(s);
    while (std::getline(scanner, nextComponent, ' ')) {
        if (nextComponent.empty()) {
            continue;
        }
        res.push_back(nextComponent);
    }
    return res;
}

void removeCharacters(std::string &s, const std::set<char> &charactersToRemove) {
    s.erase(std::remove_if(s.begin (), s.end (), [&charactersToRemove](const char &c) {
        return charactersToRemove.count(c) > 0;
    }), s.end());
}

std::set<PeakId> getDuplicates(std::deque<Peak> peaks, const std::string &mountsDatabase,
        std::set<PeakId> &savedIds)
{
    std::set<PeakId> duplicateIds;

    for (auto &peak : peaks) {
        std::transform(peak.name.begin(), peak.name.end(), peak.name.begin(),
                [](unsigned char c){ return std::tolower(c); });
    }

    std::set<char> punctuationToRemove{'(', ')', '.', '-', ','};
    auto freqWords = readStringsFromFile(mountsDatabase);
    std::set<std::string> freqSet(freqWords.begin(), freqWords.end());

    for (auto &peak : peaks) {
        auto splitName = split(peak.name);
        std::string trimmedName;
        std::deque<std::string> omittedWords;
        for (auto &word : splitName) {
            removeCharacters(word, punctuationToRemove);
            if (word.empty()) {
                continue;
            }
            if (freqSet.count(word)) {
                omittedWords.push_back(word);
            } else {
                if (trimmedName.empty()) {
                    trimmedName = word;
                } else {
                    trimmedName.append(" " + word);
                }
            }
        }
        if (trimmedName.empty()) {
            if (omittedWords.empty()) {
                trimmedName = peak.name;
            } else {
                for (const auto &w : omittedWords) {
                    trimmedName.append(w);
                    omittedWords.clear();
                }
            }
        }
        peak.name = trimmedName;
    }

    std::sort(peaks.begin(), peaks.end(), [](const Peak &f, const Peak &s) {
        return f.name < s.name;
    });

    PeakTree tree(peaks);

    struct SameNamePeaks {
        std::string name;
        std::vector<size_t> indexes;
    };

    std::deque<SameNamePeaks> duplicateNamePeaks;
    SameNamePeaks current{peaks.front().name, {0}};

    for (size_t i = 1; i < peaks.size(); ++i) {
        if (peaks[i].isExtra()) {
            continue;
        }
        const auto &name = peaks[i].name;
        if (name != current.name) {
            if (current.indexes.size() > 1) {
                duplicateNamePeaks.push_back(current);
            }
            current.name = name;
            current.indexes.clear();
        }
        current.indexes.push_back(i);
    }

    std::sort(duplicateNamePeaks.begin(), duplicateNamePeaks.end(), [](const SameNamePeaks &f, SameNamePeaks &s) {
        return f.indexes.size() > s.indexes.size();
    });

    flog() << duplicateNamePeaks.size() << " popular names\n";
    for (const auto &namePeaks : duplicateNamePeaks) {
        std::map<size_t, Geodesic> duplicates;
        for (const auto &index : namePeaks.indexes) {
            const auto &peak = peaks[index];
            duplicates[index] = Geodesic{peak.latitude, peak.longitude};
        }

        for (const auto &i : namePeaks.indexes) {
            auto it = duplicates.find(i);
            if (it == duplicates.end()) {
                continue;
            }
            auto indexGeo = *it;
            const auto &peak = peaks[i];
            std::deque<size_t> duplicatesToRemove;
            for (auto &otherIndexGeo : duplicates) {
                if (otherIndexGeo.first == i) {
                    continue;
                }
                auto dist = DistanceTools::distanceOnEarth(indexGeo.second, otherIndexGeo.second);
                if (dist > 5) {
                    continue;
                }
                const auto &otherPeak = peaks[otherIndexGeo.first];

                size_t lessProminentIndex = peak.prominence < otherPeak.prominence ? i : otherIndexGeo.first;
                size_t moreProminentIndex = lessProminentIndex == i ? otherIndexGeo.first : i;
                auto &lessProminentPeak = peaks[lessProminentIndex];
                auto &moreProminentPeak = peaks[moreProminentIndex];

                if (savedIds.count(lessProminentPeak.id)) {
                    if (savedIds.count(moreProminentPeak.id)) {
                        continue;
                    } else {
                        tree.switchNames(lessProminentIndex, moreProminentIndex);
                    }
                }

                size_t popularity = duplicates.size();

                Distance cutoffDistance;

                if (lessProminentPeak.prominence >= 100) {
                    cutoffDistance = 0.2;
                } else if (lessProminentPeak.prominence >= 50) {
                    cutoffDistance = 1;
                } else {
                    cutoffDistance = 2.5;
                }

                if (popularity < 5) {
                    cutoffDistance *= 2;
                }

                if (dist <= cutoffDistance) {
                    duplicatesToRemove.push_back(lessProminentIndex);
                    duplicateIds.insert(lessProminentPeak.id);
                    size_t moreProminentIndex = lessProminentIndex == i ? otherIndexGeo.first : i;
                    const auto &moreProminentPeak = peaks[moreProminentIndex];
                    flog() << "removing " << lessProminentPeak.id << " as double of " << moreProminentPeak.id << " both named: " << lessProminentPeak.name << " lesser prominence: " << lessProminentPeak.prominence << " dist: " << dist << " popularity: " << popularity << " cutoff: " << cutoffDistance << std::endl;
                    if (lessProminentIndex == i) {
                        break;
                    }
                }
            }
            for (const auto &indexToRemove : duplicatesToRemove) {
                duplicates.erase(indexToRemove);
            }
        }
    }
    flog() << "duplicates:" << duplicateIds.size() << std::endl;
    return duplicateIds;
}

struct RemoveDuplicatesProcessor : Processor {
    Path mountsDatabase;
    Path vipIdsDatabase;

    RemoveDuplicatesProcessor(const Path &mountsDatabase, const Path &vipIdsDatabase)
            : mountsDatabase(mountsDatabase),  vipIdsDatabase(vipIdsDatabase) {}

    std::string getName() final { return "remove_duplicates"; }

    void process(std::deque<Peak> &peaks) final {
        auto idPeaks = loadPeaks(vipIdsDatabase);
        std::set<PeakId> ids;
        for (const auto &idPeak : idPeaks) {
            ids.insert(idPeak.id);
        }
        idPeaks.clear();

        auto duplicateIds = getDuplicates(peaks, mountsDatabase, ids);

        PeakTree tree(peaks);
        std::deque<size_t> duplicateIndexes;
        for (const auto &id : duplicateIds) {
            duplicateIndexes.push_back(tree.idToIndex[id]);
        }
        std::sort(duplicateIndexes.begin(), duplicateIndexes.end(), [&peaks](size_t f, size_t s) {
            return peaks[f].definitiveAltitudeLess(peaks[s]);
        });
        assert(peaks[duplicateIndexes.front()].altitude <= peaks[duplicateIndexes.back()].altitude);

        for (size_t i = 1; i < duplicateIndexes.size(); ++i) {
            assert(duplicateIndexes[i] != duplicateIndexes[i-1]);
        }
        std::set<size_t> duplicateSet(duplicateIndexes.begin(), duplicateIndexes.end());
        assert(duplicateSet.size() == duplicateIndexes.size());
        for (size_t i : duplicateIndexes) {
            tree.unnamePeak(i);
        }
        removeZeroMarkedPeaks(peaks);
    }
};

struct LowerUnnamedProcessor : Processor {
    Altitude prominenceCutoff = 100;

    LowerUnnamedProcessor(Altitude prominenceCutoff = 100) : prominenceCutoff(prominenceCutoff) {}

    std::string getName() final { return "lower_unnamed"; }

    void process(std::deque<Peak> &peaks) final {
        PeakTree tree(peaks);
        tree.resolveSameAltitudeParenting();
        tree.lowerUnnamed(prominenceCutoff);
        removeZeroMarkedPeaks(peaks);
    }
};

struct AdjustAltitudeProcessor : Processor {
    Path fixedAltitudeIdsPath;
    Path refrencePeaksPath;

    AdjustAltitudeProcessor(Path fixedAltitudeIdsPath, Path refrencePeaksPath)
        : fixedAltitudeIdsPath(std::move(fixedAltitudeIdsPath)), refrencePeaksPath(std::move(refrencePeaksPath)) {}

    std::string getName() final { return "adjust_altitude"; }

    void process(std::deque<Peak> &peaks) final {
        auto idPeaks = loadPeaks(fixedAltitudeIdsPath);
        std::set<PeakId> fixedAltitudeIds;
        for (const auto &idPeak : idPeaks) {
            fixedAltitudeIds.insert(idPeak.id);
        }
        auto refPeaks = loadPeaks(refrencePeaksPath);
        sortById(refPeaks);

        for (auto &peak : peaks) {
            if (peak.altitude < 0) {
                peak.altitude = 0; // TODO: move to addMissing
            }
        }

        sortById(peaks);
        PeakTree tree(peaks);

        tree.checkConsistency();

        struct DiffPeak {
            PeakId id;
            std::string name;
            Geodesic geo;

            Altitude demAltitude;
            Altitude sourceAltitude;

            Altitude diff() const {
                return sourceAltitude - demAltitude;
            }
            double proportionalDiff() const {
                return static_cast<double>(sourceAltitude) / demAltitude;
            }
            bool shouldFix() const {
                if (sourceAltitude == 0) {
                    return false;
                }
                return diff() > 0 && proportionalDiff() < 2;
            }
        };
        std::deque<DiffPeak> diffPeaks;
        auto lambda = [&](size_t index, size_t refIndex, DualState state) {
            if (state == kDualStateBoth) {
                auto &peak = peaks.at(index);
                auto &refPeak = refPeaks.at(refIndex);
                auto diffPeak = DiffPeak{peak.id, peak.name, peak.geo(), peak.altitude, refPeak.altitude};
                if (diffPeak.diff() > 0) {
                    diffPeaks.emplace_back(diffPeak);
                }

                if (diffPeak.shouldFix() || (refPeak.altitude > 0 && fixedAltitudeIds.count(peak.id))) {
                    tree.changeAltitude(peak.id, refPeak.altitude);
                }
            }
        };

        dualIteration(peaks, refPeaks, lambda, true);

        std::sort(diffPeaks.begin(), diffPeaks.end(), [](const DiffPeak &lhs, const DiffPeak &rhs) {
            return lhs.proportionalDiff() > rhs.proportionalDiff();
        });

        for (const auto &diffPeak : diffPeaks) {
            flog() << diffPeak.id << "\t" << diffPeak.name << "\t" << diffPeak.sourceAltitude <<
            "\t" << diffPeak.demAltitude << "\t" << diffPeak.proportionalDiff() << std::endl;
        }
        removeZeroMarkedPeaks(peaks);
    }
};

struct NHNProcessor : Processor {
    std::string getName() final { return "nhn"; }
    struct PeakBasicInfo {
        Geodesic geo;
        Altitude altitude;
        uint64_t index;
    };

    void process(std::deque<Peak> &peaks) final {
        std::sort(peaks.begin(), peaks.end(), [](const Peak &one, const Peak &other) {
            return one.definitiveAltitudeLess(other);
        });

        std::map<PlanarIndex, std::deque<PeakBasicInfo>> peakInfos;
        std::deque<PlanarIndex> peakPlanars;

        for (uint64_t i = 0; i < peaks.size(); ++i) {
            const auto &peak = peaks[i];
            auto planar = PlanarIndex{int(std::floor(peak.latitude)), int(std::floor(peak.longitude))};

            peakInfos[planar].emplace_back(PeakBasicInfo{{peak.latitude, peak.longitude}, peak.altitude, i});
            if (i < peaks.size() - 1) {
                peakPlanars.push_back(planar);
            }
        }
        uint64_t kNoNhn = peaks.size();

        std::array<Geodesic, 2> region{};

        int count = 0;
        double nextPercent = 0.1;
        for (const auto &planar : peakPlanars) {
            const auto &peakInfo = peakInfos[planar].front();

            assert(static_cast<uint64_t>(count) == peakInfo.index);

            if ((100. * count) / peaks.size() > nextPercent) {
                flog() << nextPercent << "% complete" << std::endl;
                nextPercent += 1;
            }

            Distance minNeighborDistance = kEarthRadius * 4;
            uint64_t nhn = kNoNhn;
            double radius = 0;
            bool foundLessThanEdge = false;
            while (nhn == kNoNhn || !foundLessThanEdge) {
                Distance edgeMinDistance = kEarthRadius * 4;
                region[0] = Geodesic{planar.x - radius, planar.y - radius};
                region[1] = Geodesic{planar.x + radius + 1, planar.y + radius + 1};
                for (int i = 0; i < 2; ++i) { // const coordinate
                    for (int j = 0; j < 2; ++j) { // i value
                        for (int k = 0; k <= 2 * radius; ++k) {
                            Geodesic sideFrom{};
                            Geodesic sideTo{};
                            sideFrom.c[i] = region[j].c[i];
                            sideTo.c[i] = region[j].c[i];
                            sideFrom.c[1 - i] = region[0].c[1 - i] + k;
                            sideTo.c[1 - i] = sideFrom.c[1 - i] + 1;
                            sideFrom.normalize();
                            sideTo.normalize();
                            if (sideFrom.longitude == 180 && sideTo.longitude != 180) {
                                sideFrom.longitude = -180;
                            }
                            if (sideFrom.latitude == sideTo.latitude && std::abs(sideFrom.latitude) == 90) {
                                continue;
                            }

                            if (std::abs(sideFrom.latitude) == 90) {
                                sideFrom.longitude = sideTo.longitude;
                            }
                            if (std::abs(sideTo.latitude) == 90) {
                                sideTo.longitude = sideFrom.longitude;
                            }

                            if (sideTo.latitude <= sideFrom.latitude && sideTo.longitude <= sideFrom.longitude) {
                                std::swap(sideTo, sideFrom);
                            }

                            edgeMinDistance = std::min(DistanceTools::minMaxDistanceToQuadrant(peakInfo.geo, sideFrom, sideTo).min, edgeMinDistance);
                        }
                    }
                }

                if (minNeighborDistance < edgeMinDistance) {
                    foundLessThanEdge = true;
                }

                for (int i = -static_cast<int>(radius); i < radius || (radius == 0. && i == 0); ++i) {
                    for (int side = 0; side < 4 && (radius > 0 || side == 0); ++side) {
                        int constDeltaCoordinate = side % 2;
                        int constDeltaCoordinateSign = side / 2 ? 1 : -1;
                        PlanarIndex delta{};
                        delta.c[constDeltaCoordinate] = radius * constDeltaCoordinateSign;
                        delta.c[1 - constDeltaCoordinate] = i;
                        PlanarIndex neighborPlanar = planar + delta;
                        for (const auto &neighbor : peakInfos[neighborPlanar]) {
                            if (neighbor.index == peakInfo.index) {
                                continue;
                            }
                            if (peaks.at(neighbor.index).altitude == peaks.at(peakInfo.index).altitude) {
                                continue;
                            }
                            Distance newDistance = DistanceTools::distanceOnEarth(peakInfo.geo, neighbor.geo);
                            if (newDistance < minNeighborDistance) {
                                nhn = neighbor.index;
                                minNeighborDistance = newDistance;
                                if (minNeighborDistance < edgeMinDistance) {
                                    foundLessThanEdge = true;
                                }
                            }
                        }
                    }
                }
                ++radius;
                if (radius >= 89) {
                    nhn = kNoNhn;
                    break;
                }
            }

            if (nhn != kNoNhn) {
                auto &peak = peaks[peakInfo.index];
                const auto &nhPeak = peaks[nhn];
                peak.nhn = nhPeak.id;
                correctIsolation(peak, nhPeak);
            }
            ++count;
            peakInfos[planar].pop_front();
        }
        bruteForce(peaks);
    }

    static void bruteForce(std::deque<Peak> &peaks) {
        std::vector<size_t> notFound;
        for (auto i = 0u; i < peaks.size(); ++i) {
            if (peaks[i].nhn == kNoId) {
                notFound.emplace_back(i);
            }
        }
        for (auto toFind : notFound) {
            auto &peak = peaks[toFind];
            int64_t nhn = -1;
            Distance minDistance = 4 * kEarthRadius;
            for (int64_t i = peaks.size() - 1; i >= 0; --i) {
                if (static_cast<uint64_t>(i) == toFind) {
                    continue;
                }
                const auto &otherPeak = peaks[i];
                if (otherPeak.altitude <= peak.altitude) {
                    break;
                }

                Distance newDistance = DistanceTools::distanceOnEarth({peak.latitude, peak.longitude}, {otherPeak.latitude, otherPeak.longitude});
                if (newDistance < minDistance) {
                    nhn = i;
                    minDistance = newDistance;
                }
            }
            if (nhn == -1) {
                flog() << "cant find nhn for " << peak.toDescriptiveString() << std::endl;
            } else {
                const auto &nhPeak = peaks[nhn];

                peak.nhn = nhPeak.id;
                correctIsolation(peak, nhPeak);
            }
        }
    }
};

struct CorrectIsolationByNHNProcessor : Processor {
    std::string getName() final { return "correct_isolation_by_nhn"; }
    void process(std::deque<Peak> &peaks) final {
        PeakTree tree(peaks);
        tree.correctIsolationByNHN();
    }
};

struct FillProportionalProminenceProcessor : Processor {
    std::string getName() final { return "fill_proportial_prominence"; }
    void process(std::deque<Peak> &peaks) final {
        std::sort(peaks.begin(), peaks.end(), [](const Peak &first, const Peak &second) {
            return first.prominence > second.prominence;
        });

        PeakTree tree(peaks);

        std::function<void(Peak &)> fillPP = [&](Peak &peak) {
            if (peak.prominence == 0) {
                peak.proportionalProminence = 0;
                return;
            }
            if (peak.islandParent == kNoId) {
                return;
            }
            if (peak.proportionalProminence != kNoAltitude) {
                return;
            }
            if (peak.islandParent == peak.id) {
                assert(peak.altitude == peak.prominence);
                peak.proportionalProminence = peak.prominence;
            } else {
                if (tree.idToIndex.count(peak.islandParent) == 0) {
                    std::cout << peak.toDescriptiveString() << std::endl;
                    assert(false);
                }
                auto &parent = peaks[tree.idToIndex.at(peak.islandParent)];
                if (parent.proportionalProminence == kNoAltitude) {
                    fillPP(parent);
                }
                if (parent.proportionalProminence != kNoAltitude) {
                    const auto keycolAlt = peak.altitude - peak.prominence;
                    peak.proportionalProminence = Altitude(double(parent.proportionalProminence) / (parent.altitude - keycolAlt) * peak.prominence);
                }
            }
        };

        for (auto &peak : peaks) {
            fillPP(peak);
        }
    }
};

struct AddPassesProcessor : Processor {
    Path passesPath;

    explicit AddPassesProcessor(Path passesPath) : passesPath(std::move(passesPath)) {}

    std::string getName() final { return "add_passes"; }
    void process(std::deque<Peak> &peaks) final {
        addPasses(peaks, passesPath);
    }
};

std::vector<std::unique_ptr<Processor>> getPostProcessors(const Path &resourcesDir, const Path &rawPeaksPath) {
    ResourcesPaths res{resourcesDir};
    std::vector<std::unique_ptr<Processor>> result;
    result.emplace_back(std::make_unique<RemoveDuplicateIdsProcessor>());
    result.emplace_back(std::make_unique<AddMissingPeaksProcessor>(rawPeaksPath));
    result.emplace_back(std::make_unique<RemoveDuplicatesProcessor>(res.mounts(), res.vipIds()));
    result.emplace_back(std::make_unique<LowerUnnamedProcessor>());
    result.emplace_back(std::make_unique<AdjustAltitudeProcessor>(res.fixedAltitude(), rawPeaksPath));
    result.emplace_back(std::make_unique<LowerUnnamedProcessor>());
    result.emplace_back(std::make_unique<NHNProcessor>());
    result.emplace_back(std::make_unique<CorrectIsolationByNHNProcessor>());
    result.emplace_back(std::make_unique<FillProportionalProminenceProcessor>());
    result.emplace_back(std::make_unique<AddPassesProcessor>(res.passes()));
    return result;
}

struct Pipeline {
    static constexpr auto kInputPeaksFilename = "input_peaks.txt";
    static constexpr auto kRefinedPeaksFilename = "refined_peaks.txt";
    static constexpr auto kPeaksWithProminenceFileneme = "peaks_with_prominence.txt";

    std::optional<Path> generationDir;
    std::optional<Path> referenceDir;

    std::vector<Path> generatedFilenames;
    std::vector<std::unique_ptr<Processor>> processors;

    explicit Pipeline(std::optional<Path> generationDir = std::nullopt, std::optional<Path> referenceDir = std::nullopt)
        : generationDir(std::move(generationDir)), referenceDir(std::move(referenceDir)) {}

    void applyProcessors(std::deque<Peak> &peaks, bool overwrite = true) {
        for (auto &processor : processors) {
            flog() << "applying " << processor->getName() << std::endl;

            std::stringstream nextFilenameStream;
            nextFilenameStream << std::hex << generatedFilenames.size() << "_" << processor->getName() + ".txt";
            auto nextFilename = generationDir ? std::optional<Path>(nextFilenameStream.str()) : std::nullopt;
            auto nextPath = nextFilename ? std::optional<Path>(*generationDir / *nextFilename) : std::nullopt;

            if (referenceDir && !generatedFilenames.empty()) {
                auto refPath = *referenceDir / generatedFilenames.back();
                peaks = loadPeaks(refPath);
            }

            if (!overwrite && nextPath && std::filesystem::exists(*nextPath)) {
                flog() << "already have " << *nextPath << std::endl;
                peaks = loadPeaks(*nextPath);
                generatedFilenames.emplace_back(*nextFilename);
                continue;
            }

            try {
                processor->process(peaks);
            } catch (const std::logic_error &error) {
                std::cout << "err: " << error.what() << std::endl;
                if (generationDir) {
                    savePeaks(peaks, *generationDir / "with_error.txt");
                }
                throw;
            }
            if (nextPath) {
                savePeaks(peaks, *nextPath);
                generatedFilenames.emplace_back(*nextFilename);
            }
        }
    }

    void addPostProcessing(const Path &resourcesDir, const Path &rawInputPeaksPath) {
        auto postProcessors = getPostProcessors(resourcesDir, rawInputPeaksPath);
        processors.insert(processors.end(), std::make_move_iterator(postProcessors.begin()), std::make_move_iterator(postProcessors.end()));
    }
};

inline void postProcess(std::deque<Peak> &peaks, const Path &resourcesDir, Path rawInputPeaksPath,
                        std::optional<Path> generationDir = std::nullopt, std::optional<Path> referenceDir = std::nullopt, bool overwrite = true)
{
    Pipeline pipeline{generationDir, referenceDir};
    pipeline.addPostProcessing(resourcesDir, rawInputPeaksPath);
    pipeline.applyProcessors(peaks, overwrite);
}

inline void postProcess(Path inputPeaksFile, Path outputPeaksFile, Path resourcesDir, Path rawInputPeaksPath,
                        std::optional<Path> generationDir = std::nullopt, std::optional<Path> referenceDir = std::nullopt, bool overwrite = true)
{
    auto peaks = loadPeaks(inputPeaksFile);
    postProcess(peaks, resourcesDir, rawInputPeaksPath, generationDir, referenceDir, overwrite);
    savePeaks(peaks, outputPeaksFile);
}

#endif // PIPELINE_H
