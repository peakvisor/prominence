#ifndef DEM_PIPELINE_H
#define DEM_PIPELINE_H

#include "pipeline.h"
#include "prominence_calculation.h"

template <typename FileFormat>
struct ProminenceCalculationProcessor : Processor {
    Path demsDirectory;
    Rectangle region;
    Path backUpPeaksFile;

    ProminenceCalculationProcessor(Path demsDirectory, const Rectangle &region, Path backUpPeaksFile)
    : demsDirectory(std::move(demsDirectory)), region(region), backUpPeaksFile(std::move(backUpPeaksFile)) {}

    std::string getName() final { return "prominence_calculation"; }

    void process(std::deque<Peak> &peaks) final {
        auto dem = FileFormat::getDEM(demsDirectory, region);
        std::deque<Peak> filteredByDem;
        for (auto &peak : peaks) {
            if (dem.geodesicIsInside(peak.geo())) {
                filteredByDem.emplace_back(peak);
            }
        }
        peaks = calculateProminences(dem, filteredByDem, backUpPeaksFile);
    }
};

template <typename FileFormat>
void calculateProminences(Path demsDirectory, const Rectangle &region, Path inputPeaksFile, Path outputPeaksFile) {
    std::cout << "calc prom\n";
    auto processor = ProminenceCalculationProcessor<FileFormat>(demsDirectory, region, outputPeaksFile);
    auto peaks = loadPeaks(inputPeaksFile);
    processor.process(peaks);
    std::sort(peaks.begin(), peaks.end(), [](const Peak &lhs, const Peak &rhs) {
        return lhs.prominence > rhs.prominence;
    });
    savePeaks(peaks, outputPeaksFile);
}

#endif // DEM_PIPELINE_H
