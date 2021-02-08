#define DEBUG_ASSERTS

#include "distance.h"
#include "distance_tests.h"
#include "../peaks/peaks_processing.h"
#include "../pipeline/pipeline.h"
#include "prominence_calculation.h"
#include "geomorphology.h"
#include "../pipeline/dem_pipeline.h"
#include "../peaks_on_dem_processing/peaks_on_dem_processing.h"
#include "../dem/dem_processing.h"

#include <random>
#include <chrono>

std::vector<Rectangle> splitDemsIntoRegions(const std::filesystem::path &dir) {
    std::vector<Rectangle> result;
    auto joinRectangles = [&result](auto lIndex, auto rIndex) {
        auto joined = rectangularShell(result.at(lIndex), result.at(rIndex));
        result.at(std::min(lIndex, rIndex)) = joined;
        result.at(std::max(lIndex, rIndex)) = result.back();
        result.pop_back();
    };

    auto tryJoin = [&result, &joinRectangles](auto &lIndex) {
        for (auto rIndex = 0u; rIndex < result.size(); ++rIndex) {
            if (rIndex == lIndex) {
                continue;
            }
            if (rectanglesNear(result.at(lIndex), result.at(rIndex))) {
                joinRectangles(lIndex, rIndex);
                return true;
            }
        }
        return false;
    };

    for (auto &dirIt : std::filesystem::directory_iterator(dir)) {
        auto filename = dirIt.path().filename();
        PlanarIndex planarIndex = DEMFilenameFormatter::parseLatLng(filename);
        Rectangle rectangle{planarIndex, planarIndex + kDiagonalShift};
        result.emplace_back(rectangle);
        bool joinedSomething;
        do {
            joinedSomething = false;
            for (auto i = 0u; i < result.size(); ++i) {
                if (tryJoin(i)) {
                    joinedSomething = true;
                    break;
                }
            }
        } while (joinedSomething);
    }
    return result;
}

template <typename FileFormat>
void testRegion(const Path &pipelineTestingDir, bool overwriteReference) {
    auto resourcesDir = pipelineTestingDir / "res";
    if (!std::filesystem::exists(resourcesDir)) {
        std::cout << "no resources dir found at " << resourcesDir << std::endl;
        return;
    }
    auto demsDir = pipelineTestingDir / (FileFormat::prefix() + "dems");
    if (!std::filesystem::exists(demsDir)) {
        std::cout << demsDir << "no dems dir found at " << demsDir << std::endl;
        return;
    }

    for (auto &region : splitDemsIntoRegions(demsDir)) {
        auto regionName = FileFormat::prefix() +
            (region.end - region.origin == kDiagonalShift
            ? DEMFilenameFormatter::formatLatLng(region.origin)
            : DEMFilenameFormatter::formatLatLng(region.origin) + "-" + DEMFilenameFormatter::formatLatLng(region.end - kDiagonalShift));
        std::cout << "running tests for " << regionName << std::endl;

        auto regionDir = pipelineTestingDir / regionName;
        auto referenceDir = regionDir / "ref";
        bool isTesting = std::filesystem::exists(referenceDir) && !overwriteReference;
        auto generationDir = isTesting ? regionDir / "test" : referenceDir;
        if (!std::filesystem::exists(generationDir)) {
            std::filesystem::create_directories(generationDir);
        }
        if (!isTesting) {
            filterByCoordinates(resourcesDir / Pipeline::kRefinedPeaksFilename, referenceDir / Pipeline::kRefinedPeaksFilename,
                {region.origin.x - 0.05, region.origin.y - 0.05}, {region.end.x + 0.05, region.end.y + 0.05});
            filterByCoordinates(resourcesDir / Pipeline::kInputPeaksFilename, referenceDir / Pipeline::kInputPeaksFilename,
                {region.origin.x - 0.05, region.origin.y - 0.05}, {region.end.x + 0.05, region.end.y + 0.05});
        }

        Pipeline pipeline(generationDir, isTesting ? std::optional<Path>(referenceDir) : std::nullopt);
        Path backupProminenceCalculationPath = generationDir /  Pipeline::kPeaksWithProminenceFileneme;
        pipeline.processors.emplace_back(std::make_unique<ProminenceCalculationProcessor<FileFormat>>(demsDir, region, backupProminenceCalculationPath));
        pipeline.addPostProcessing(resourcesDir, referenceDir / Pipeline::kInputPeaksFilename);
        auto processorsCount = pipeline.processors.size();
        auto peaks = loadPeaks(referenceDir / Pipeline::kRefinedPeaksFilename);
        pipeline.applyProcessors(peaks);

        if (isTesting) {
            size_t totalChecked = 0;
            for (const auto &peaksFilename : pipeline.generatedFilenames) {
                auto referencePath = referenceDir / peaksFilename;
                auto generatedPath = generationDir / peaksFilename;
                if (std::filesystem::exists(referencePath) && std::filesystem::exists(generatedPath)) {
                    if (!peakFilesAreTheSame(referencePath, generatedPath, true)) {
                        throw std::logic_error("peaks files differ: " + std::string(referencePath) + " and " + std::string(generatedPath));
                    }
                    ++totalChecked;
                }
            }
            if (totalChecked != processorsCount) {
                throw std::logic_error("some files missing");
            }
        }
        assert(peakFilesAreTheSame(backupProminenceCalculationPath, generationDir / pipeline.generatedFilenames.front()));
    }
}

void testPipeline(const std::filesystem::path &testDir, bool overwriteReference = false) {
    testRegion<SDFileFormat>(testDir, overwriteReference);
    testRegion<HDFileFormat>(testDir, overwriteReference);
}

void testPostProcessingResult(const Path &resourcesDir, const Path &tmpDir,
        const Path &srcPeaksPath, const Path &postProcessedPath)
{
    Pipeline pipeline{tmpDir};
    pipeline.addPostProcessing(resourcesDir, resourcesDir / "input_peaks.txt");
    auto peaks = loadPeaks(srcPeaksPath);
    pipeline.applyProcessors(peaks);
    std::cout << "peaks files are the same: " << peakFilesAreTheSame(postProcessedPath, tmpDir / pipeline.generatedFilenames.back()) << std::endl;
}

template <typename FileFormat>
void getTilesWithoutBorder(const Path &demsDir) {
    std::vector<std::pair<int, Path>> zeroBound;

    int index = 0;
    for (auto &dirIt : std::filesystem::directory_iterator(demsDir)) {
        auto filename = dirIt.path().filename();
        if (!filename.has_extension()) {
            continue;
        }
        ++index;
        PlanarIndex p = DEMFilenameFormatter::parseLatLng(filename);
        auto dem = FileFormat::getDEM(demsDir, {p, p + kDiagonalShift});
        bool hasZeroBorder = true;
        for (int i : {0, 1}) {
            for (Index x : {0, FileFormat::cellsResolution.c[i] - 1}) {
                for (Index y = 0; y < FileFormat::cellsResolution.c[1 - i]; ++y) {
                    auto bp = planarFromIndexedCoordinates(i, x, y);
                    if (dem.get(bp) != kOceanFloor)
                    {
                        hasZeroBorder = false;
                        break;
                    }
                }
                if (!hasZeroBorder) {
                    break;
                }
            }
            if (!hasZeroBorder) {
                break;
            }
        }
        if (!hasZeroBorder) {
            continue;
        }

        int pozitiveCount = 0;
        dem.doForEachPlanar([&dem, &pozitiveCount](const PlanarIndex &p) {
            if (dem.get(p) != kOceanFloor) { ++pozitiveCount; }
        });
        zeroBound.emplace_back(pozitiveCount, filename);
    }
    std::sort(zeroBound.rbegin(), zeroBound.rend());
    for (int i = 0; i < 10; ++i) {
        std::cout << zeroBound[i].second << " " << zeroBound[i].first << std::endl;
    }
}

int main(int argc, char * argv[]) {
    if (argc == 1234) {
        applyProminenceCalculation(argc, argv);
    } else if (argc == 12345) {
        applyPeakListTools(argc, argv);
    } else if (argc == 123456) {
        applyGeomorphologyTools(argc, argv);
    } else if (argc == 343) {
        applyPeaksOnDemProcessing(argc, argv);
    } else if (argc == 432) {
        applyDEMProcessing(argc, argv);
    }
    const auto t1 = std::chrono::high_resolution_clock::now();
    if (argc <= 1) {
        std::cout << "no test dir specified" << std::endl;
        return 0;
    }

    Path testDir{argv[1]};
    FileLog::setUpFileLog(testDir / "log.txt");
    testDistance();
    testPipeline(testDir);
//    getTilesWithoutBorder<SDFileFormat>("/Users/ka/pc/prod/sd_dems");
//    calcRidges({33, 126}, {33, 126}, "/Users/ka/pc/prod/hd_dems",
//            "/Users/ka/pc/prod/test/hd_N33E126/ref/0_prominence_calculation.txt",
//            normalizedDir("/Users/ka/pc/prod/test/ridges"), 0);

    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed << "finished in " << ms.count() / 1000. << "s\n";
}