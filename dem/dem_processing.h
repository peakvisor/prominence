#ifndef DEM_PROCESSING_H
#define DEM_PROCESSING_H

#include "common.h"
#include "index.h"
#include "field_view.h"
#include "dem_filename_formatter.h"
#include "dem_grid.h"
#include "logging.hpp"

template <typename ValueType>
void serializeNonEmptyVector(const std::vector<ValueType> &container, const std::string &filePath) {
    std::ofstream output{filePath};
    if (!output) {
        std::cout << "\tfailed to serialize to '" << filePath << "'!\n";
        assert(false);
    }
    if (container.empty()) {
        std::cout << "\tattempting to serialize empty container to '" << filePath << "'!\n";
        assert(false);
        return;
    }
    output.write(reinterpret_cast<const char *>(&container.front()), container.size() * sizeof(container[0]));
}

static const int kTile = 2, kTileArea = kTile * kTile;
void prepareFieldForCompression(const FieldView &input, std::vector<Altitude> &output, std::vector<Altitude> &auxiliary) {
    int downscaledLength = static_cast<int>(input.area() / kTileArea);

    input.tileInto(auxiliary, kTile);
    for (int i = 0, iin = 0, iout = downscaledLength; i < downscaledLength; ++i, iin += kTileArea, iout += kTileArea - 1) {
        auto reference = output[i] = auxiliary[iin];
        for (int j = 1; j < kTileArea; ++j) {
            output[iout + j - 1] = auxiliary[iin + j] - reference;
        }
    }
}

std::string preparedFileName(const PlanarIndex &latLng, const std::string &prefix = "tHDDEM_c") {
    static const Index kCurrentVersion = 0;
    static const std::string kCurrentFormat = "M2SEP";
    return prefix + DEMFilenameFormatter::formatLatLng(latLng) + "_v" + std::to_string(kCurrentVersion) + "_f" + kCurrentFormat + ".edem";
}

void prepareForCompressionNDEMFiles(const std::string &inputDirectory,
    const std::string &outputDirectory, const PlanarIndex &from = kLatLngFrom,
    const PlanarIndex &to = kLatLngTo, int margin = 2)
{
    RegularDEMGrid grid{inputDirectory, from - PlanarIndex{1, 1}, to + PlanarIndex{1, 1}, {1200, 1200}, {1200, 1200}};
    auto extendedResolution = grid.cellsResolution.x + 2 * margin;
    std::vector<Altitude> buffer, auxiliary;
    FieldView field{buffer, extendedResolution};

    auxiliary.resize(extendedResolution * extendedResolution);

    int d__nonemptyCellsCounter = 0;
    for (Index latitude = from.latitude; latitude <= to.latitude; ++latitude) {
        for (Index longitude = from.longitude; longitude <= to.longitude; ++longitude) {
            PlanarIndex latLng{latitude, longitude};
            ++d__nonemptyCellsCounter;

            auto origin = grid.geodesicToPlanar({double(latitude), double(longitude)});
            origin.x -= margin; origin.y -= margin;

            bool d__nontrivial = false;
            for (Index x = 0; x < extendedResolution; ++x) {
                for (Index y = 0; y < extendedResolution; ++y) {
                    auto h = field.get({x, y}) = grid.get({origin.x + x, origin.y + y});
                    d__nontrivial = d__nontrivial || h; // we suppose kOceanLevel/... is 0
                }
            }
            if (!d__nontrivial) {
                continue;
            }

            prepareFieldForCompression(field, buffer/* 'inplace' */, auxiliary);
            serializeNonEmptyVector(buffer, outputDirectory + preparedFileName(latLng));
        }
    }
}

void unprepareFieldForCompression(const std::vector<Altitude> &input, FieldView &output, std::vector<Altitude> &auxiliary) {
    auto downscaledLength = static_cast<int>(output.area() / kTileArea);

    for (int i = 0, iout = 0, iin = downscaledLength; i < downscaledLength; ++i, iout += kTileArea, iin += kTileArea - 1) {
        auto reference = auxiliary[iout] = input[i];
        for (int j = 1; j < kTileArea; ++j) {
            auxiliary[iout + j] = input[iin + j - 1] + reference;
        }
    }
    output.untileFrom(auxiliary, kTile);
}

void deEDEMFiles(const std::string &inputDirectory, const std::string &outputDirectory,
    const PlanarIndex &from = kLatLngFrom, const PlanarIndex &to = kLatLngTo)
{
    constexpr int kOriginalResolution = 1200, kMargin = 2,
        kExtendedResolution = kOriginalResolution + 2 * kMargin;
    std::vector<Altitude> buffer;
    std::vector<Altitude> auxiliary(kExtendedResolution * kExtendedResolution);
    std::vector<Altitude> supposedOriginal;
    FieldView field{buffer, kExtendedResolution};
    FieldView supposedOriginalField{supposedOriginal, kOriginalResolution};

    for (Index latitude = from.latitude; latitude <= to.latitude; ++latitude) {
        for (auto longitude = from.longitude; longitude <= to.longitude; ++longitude) {
            PlanarIndex latLng{latitude, longitude};
            auto filePath = inputDirectory + preparedFileName(latLng, "tDEM_c");
            std::cout << "trying " << filePath << std::endl;
            std::ifstream input{filePath};
            if (!input) {
                continue;
            }
            std::cout << "read: " << filePath << std::endl;
            input.read(reinterpret_cast<char *>(&buffer.front()), field.area() * sizeof(Altitude));

            unprepareFieldForCompression(buffer, field, auxiliary);
            field.copy(supposedOriginalField, kOriginalResolution, kOriginalResolution, kMargin, kMargin);
            serializeNonEmptyVector(supposedOriginal, outputDirectory + DEMFilenameFormatter::normalizedFileName(latLng));
        }
    }
}

void copyNDEMFiles(const std::string &inputDirectory, const std::string &outputDirectory,
    const PlanarIndex &from, const PlanarIndex &to)
{
    static const int kLength = 1200 * 1200;
    std::vector<Altitude> buffer(kLength);
    for (Index x = from.x; x <= to.x; ++x) {
        for (Index y = from.y; y <= to.y; ++y) {
            PlanarIndex latLng{x, y};
            std::ifstream input{inputDirectory + DEMFilenameFormatter::normalizedFileName(latLng)};
            if (!input) { continue; }
            input.read(reinterpret_cast<char *>(&buffer.front()), kLength * sizeof(Altitude));
            std::ofstream output{outputDirectory + DEMFilenameFormatter::normalizedFileName(latLng)};
            output.write(reinterpret_cast<char *>(&buffer.front()), kLength * sizeof(Altitude));
        }
    }
}

void mergeHDDEMs(const std::string &inputDirectory, const std::string &outputDirectory, const PlanarIndex &at) {
    RegularDEMGrid grid{inputDirectory, at * 3, at * 3 + PlanarIndex{2, 2}, {1200, 1200}, {1200, 1200}};
    if (grid.empty()) {
        return;
    }
    auto origin = grid.indexFrom();
    auto size = PlanarIndex{3600, 3600}; // grid includes both from and to
//    printf("merging [%d,%d]->(%d,%d) grid.nonempty=%d size=%d,%d\n", from.x, from.y, to.x, to.y, (int)grid.views.size(), size.x, size.y);

    std::vector<Altitude> buffer;
    FieldView field{buffer, size.x, size.y};
    for (Index x = 0; x < size.x; ++x) {
        for (Index y = 0; y < size.y; ++y) {
            field.get({x, y}) = grid.get({origin.x + x, origin.y + y});
        }
    }
    auto outputFilePath = outputDirectory + DEMFilenameFormatter::normalizedFileName(at);
    serializeNonEmptyVector(buffer, outputFilePath);

}

void tileNDEMs(const std::string &inputDirectory, const std::string &outputDirectory, const PlanarIndex &from = kLatLngFrom, const PlanarIndex &to = kLatLngTo, int originalFileSide = 1200, int xTileSizeLog = 5, int yTileSizeLog = 6) {
    int xTileSize = 1 << xTileSizeLog;
    int yTileSize = 1 << yTileSizeLog;

    RegularDEMGrid grid{inputDirectory, from, to, {originalFileSide, originalFileSide}, originalFileSide};
    std::vector<Altitude> buffer;

    TiledFieldView tiledView{buffer, kBigTiledFileSizeX, kBigTiledFileSizeY};

    for (int x = 0; x <= to.x - from.x; ++x) {
        for (int y = 0; y <= to.y - from.y; ++y) {
//            int viewIndex = grid.viewIndexFromCellIndex({x, y});
            auto &view = grid.views[grid.viewIndexFromCellIndex({x, y})];
            if (view.empty()) {
                continue;
            }
            view.tileInto(buffer, xTileSize, yTileSize, kBigTiledFileSizeX, kBigTiledFileSizeY);
            for (int i = 0; i < 1200; ++i) {
                for (int j = 0; j < 1200; ++j) {
                    if (view.get({i, j}) != tiledView.get({i, j})) {
                        std::cout << i << " " << j << ": " << view.get({i, j}) << " vs " << tiledView.get({i, j}) << std::endl;
                        assert(false);
                    }
                }
            }
            auto outputName = outputDirectory + DEMFilenameFormatter::normalizedFileName({x + from.x, y + from.y});
            serializeNonEmptyVector(buffer, outputName);
            std::cout << "Have tiled to: " << outputName << std::endl;
        }
    }
}

void testTiling() {
    std::string input = "/Users/ka/xcproj/prominence_calc/all/";
    std::string output = "/Users/ka/xcproj/prominence_calc/TILED/";

    PlanarIndex outputSize = {1216, 1216};
    RegularDEMGrid grid{input, {45, 8}, {45, 8},  {1200, 1200}, {1200, 1200}};
    std::vector<Altitude> buffer;
    FieldView_<AssymetricTiledIndexer<5, 6>> tiledView{buffer, outputSize.x, outputSize.y};
    auto &view = grid.views[0];
    view.tileInto(buffer, 32, 64, outputSize.x, outputSize.y);
    for (int i = 0; i < 1200; ++i) {
        for (int j = 0; j < 1200; ++j) {
            if (view.get({i, j}) != tiledView.get({i, j})) {
                std::cout << i << " " << j << ": " << view.get({i, j}) << " vs " << tiledView.get({i, j}) << std::endl;
                assert(false);
            }
        }
    }
    auto outputName = output + DEMFilenameFormatter::normalizedFileName({45, 8});
    std::cout << "buffer size: " << buffer.size() << std::endl;
    serializeNonEmptyVector(buffer, outputName);

    TiledDEMGrid tiledRegion{output, {45, 8}, {45, 8}, outputSize, 1200};
    assert(grid.views.size() == 1 && tiledRegion.views.size() == 1);

    auto &sameTiledView = tiledRegion.views[0];

    for (int i = 0; i < 1200; ++i) {
        for (int j = 0; j < 1200; ++j) {
            if (view.get({i, j}) != sameTiledView.get({i, j})) {
                std::cout << i << " " << j << ": " << view.get({i, j}) << " vs " << sameTiledView.get({i, j}) << std::endl;
                assert(false);
            }
        }
    }
    std::cout << "tested tiling\n";
}

void testLoad() {
    for (int i = 0; i < 100; ++i) {
        std::cout << i << std::endl;
        MappedFile file{"/Users/ka/pc/prod/all/N28E089.ndem"};
        if (file.empty()) {
            std::cout << "empty\n";
            return;
        }
        auto heightfield = static_cast<Altitude *>(const_cast<void *>(file.contents()));
        std::vector<Altitude> v;
        for (int index = 0; index < 1200 * 1200; ++index) {
            v.emplace_back(heightfield[index]);
        }
        std::cout << v[v.size() - 1] << std::endl;
    }
}

void deEdemMergeTile(const std::string &inputDirectory, const std::string &outputDirectory, const PlanarIndex &at) {
    auto name = DEMFilenameFormatter::normalizedFileName(at);
    std::cout << "processing " << name << std::endl;
    PlanarIndex from = at * 3;
    PlanarIndex to = at * 3 + PlanarIndex{2, 2};
    static const int kOriginalResolution = 1200, kMargin = 2, kExtendedResolution = kOriginalResolution + 2 * kMargin;
    std::vector<Altitude> buffer;
    std::vector<Altitude> auxiliary(kExtendedResolution * kExtendedResolution);
    std::vector<Altitude> supposedOriginal;
    FieldView field{buffer, kExtendedResolution};
    FieldView supposedOriginalField{supposedOriginal, kOriginalResolution};

    std::vector<Altitude> bigBuffer;
    TiledFieldView bigField{bigBuffer, 3616, 3648};
    bool nonZero = false;
    for (Index latitude = from.latitude; latitude <= to.latitude; ++latitude) {
        for (auto longitude = from.longitude; longitude <= to.longitude; ++longitude) {
            PlanarIndex latLng{latitude, longitude};
            auto filePath = DEMFilenameFormatter::lastVersionEdemPath(inputDirectory,
                DEMFilenameFormatter::formatLatLng(latLng));
            std::cout << "trying " << filePath << std::endl;
            std::ifstream input{filePath};
            if (!input) {
                continue;
            }
            std::cout << "read: " << filePath << std::endl;
            input.read(reinterpret_cast<char *>(&buffer.front()), field.area() * sizeof(Altitude));

            unprepareFieldForCompression(buffer, field, auxiliary);
            field.copy(supposedOriginalField, kOriginalResolution, kOriginalResolution, kMargin, kMargin);
            int shiftLat = latitude - from.latitude;
            int shiftLong = longitude - from.longitude;

            for (Index x = 0; x < kOriginalResolution; ++x) {
                for (Index y = 0; y < kOriginalResolution; ++y) {
                    auto value = supposedOriginalField.get({x, y});
                    if (!nonZero && value != 0) {
                        nonZero = true;
                    }
                    bigField.get({x + shiftLat * kOriginalResolution, y + shiftLong * kOriginalResolution}) = value;
                }
            }
        }
    }

    if (nonZero) {
        std::cout << name << " ready\n";
        serializeNonEmptyVector(bigBuffer, outputDirectory + DEMFilenameFormatter::normalizedFileName(at));
    } else {
        std::cout << name << " skipped\n";
    }
}

bool deTile(const std::string &inputDirectory, const std::string &outputDirectory, const PlanarIndex &at) {
    Index size = 3600;
    TiledDEMGrid tiledRegion{inputDirectory, at, at, {kBigTiledFileSizeX, kBigTiledFileSizeY}, size};
    if (tiledRegion.filesCount == 0) {
        return false;
    }
    std::vector<Altitude> buffer;

    FieldView view(buffer, size);

    for (Index la = 0; la < size; ++la) {
        for (Index lo = 0; lo < size; ++lo) {
            PlanarIndex p{la, lo};
            view.get(p) = tiledRegion.get(p);
        }
    }
    auto filename = DEMFilenameFormatter::normalizedFileName(at);
    auto path = outputDirectory + filename;
    std::cout << path << std::endl;
    serializeNonEmptyVector(buffer, path);
    return true;
}

bool tile(const std::string &inputDirectory, const std::string &outputDirectory, const PlanarIndex &at) {
    Index size = 3600;
    RegularDEMGrid region{inputDirectory, at, at, {3600, 3600}, 3600};
    if (region.filesCount == 0) {
        return false;
    }
    std::vector<Altitude> buffer;

    TiledFieldView view(buffer, kBigTiledFileSizeX, kBigTiledFileSizeY);

    for (Index la = 0; la < size; ++la) {
        for (Index lo = 0; lo < size; ++lo) {
            PlanarIndex p{la, lo};
            view.get(p) = region.get(p);
        }
    }
    auto filename = DEMFilenameFormatter::normalizedFileName(at);
    serializeNonEmptyVector(buffer, outputDirectory + filename);
    return true;
}


bool insideSquare(const PlanarIndex &planar, int resolution) {
    return planar.x < resolution && planar.x >= 0 && planar.y < resolution && planar.y >= 0;
}

Geodesic geodesicFromPlanar(const PlanarIndex &subIndex, const PlanarIndex &superIndex, double invRes) {
    return {superIndex.x + subIndex.x * invRes, superIndex.y + subIndex.y * invRes};
}

void findNans(const std::string &directory, std::ostream& out, int resolution, Altitude diffCutoff) {
    std::vector<Altitude> buffer;
    TiledFieldView field{buffer, 3616, 3648};

    double invRes = 1. / resolution;

    for (int la = -90; la <= 90; ++la) {
        for (int lo = -180; lo <= 180; ++lo) {
            std::string filename = DEMFilenameFormatter::normalizedFileName(PlanarIndex{la, lo});
            std::ifstream input{directory + filename};
            if (!input) {
                continue;
            }
            std::cout << "read: " << filename << std::endl;
            input.read(reinterpret_cast<char *>(&buffer.front()), field.area() * sizeof(Altitude));

            int low = 0;
            int high = 0;
            int bigDiff = 0;

            for (int i = 0; i < resolution; ++i) {
                for (int j = 0; j < resolution; ++j) {
                    PlanarIndex p{i, j};
                    Altitude a = field.get(p);
                    Geodesic geo = geodesicFromPlanar(p, {la, lo}, invRes);
                    if (a > 10000) {
                        ++high;
                        out << "high " << geo << " a: " << a << std::endl;
                    }
                    if (a < -100) {
                        ++low;
                        out << "low  " << geo << " a: " << a << std::endl;
                    }
                    int diff = 0;
                    for (const auto &delta : kNeighborOffsets) {
                        auto neighbor = p + delta;
                        if (insideSquare(neighbor, resolution)) {
                            int nextDiff = a - field.get(neighbor);
                            if (nextDiff > diff) {
                                diff = nextDiff;
                            }
                        }
                    }
                    if (diff > diffCutoff) {
                        ++bigDiff;
                        out << "diff " << geo << " d: " << diff << std::endl;
                    }
                }
            }
            if (low > 0 || high > 0 || bigDiff > 0) {
                out << "file: " << filename << " high: " << high << " low: " << low << " diff: " << bigDiff << std::endl;
            }

        }
    }

}

void uncompress(const std::string &inputPath, const std::string &outputPath) {
    int resolution = 1200;
    int kExtendedResolution = resolution + 4;
    std::vector<Altitude> buffer;
    std::vector<Altitude> auxiliary;//(kExtendedResolution * kExtendedResolution);
    std::vector<Altitude> supposedOriginal;
    FieldView supposedOriginalField{supposedOriginal, resolution};
    FieldView field{buffer, kExtendedResolution};
    FieldView auxField{auxiliary, kExtendedResolution};
    std::ifstream input{inputPath};
    if (!input) {
        return;
    }
    std::cout << "read: " << inputPath << std::endl;
    input.read(reinterpret_cast<char *>(&buffer.front()), field.area() * sizeof(Altitude));

    unprepareFieldForCompression(buffer, field, auxiliary);
    field.copy(supposedOriginalField, resolution, resolution, 2, 2);
    serializeNonEmptyVector(supposedOriginal, outputPath);
}

void sliceAndCompress(const std::string &inputDir, const std::string &outputDir, const PlanarIndex &at) {
    Index size = 3600;
    Index smallSize = size / 3;
    TiledDEMGrid tiledRegion{inputDir, at - kDiagonalShift, at + kDiagonalShift, {kBigTiledFileSizeX, kBigTiledFileSizeY}, size};

    Index kMargin = 2;
    Index sliceSide = smallSize + kMargin * 2;

    std::vector<Altitude> buffer, aux;
    FieldView view(buffer, sliceSide);
    PlanarIndex sliceOrigin;
    bool isEmptyFlag = true;
    std::function<void(const PlanarIndex&)> lambda = [&](const PlanarIndex &planar) {
        auto alt = tiledRegion.get(planar);
        if (alt != 0) {
            isEmptyFlag = false;
        }
        view.get(planar - sliceOrigin) = alt;
    };

    for (Index la = -1; la <= 3; ++la) {
        for (Index lo = -1; lo <= 3; ++lo) {
            sliceOrigin = {size + la * smallSize - kMargin, size + lo * smallSize - kMargin};
            Rectangle slice{sliceOrigin, sliceOrigin + kDiagonalShift * sliceSide};
            isEmptyFlag = true;
            slice.applyToAll(lambda);
            if (isEmptyFlag) {
                std::cout << tiledRegion.planarToGeodesic(sliceOrigin) << " is empty\n";
                continue;
            }
            prepareFieldForCompression(view, buffer, aux);
            PlanarIndex sliceCoordinates = at * 3 + PlanarIndex{la, lo};
            serializeNonEmptyVector(buffer, outputDir + preparedFileName(sliceCoordinates));
        }
    }
}

void checkSlicing(const std::string &inputDir, const std::string &outputDir, PlanarIndex at) {
    sliceAndCompress(inputDir, outputDir, at);
    deEdemMergeTile(outputDir, outputDir, at);
}

void checkTiling(const std::string &inputDir, const std::string &outputDir, PlanarIndex at) {
    deTile(inputDir, outputDir, at);
    tile(outputDir, outputDir, at);
}

void copyCorrected() {
    std::string input = "/Users/ka/pc/prod/hd/";
    std::string output = "/Users/ka/pc/prod/hd_untiled/";
    PlanarIndex from = {28, 84};
    PlanarIndex to = from;

    std::string cinput = "/Users/ka/pc/prod/hd_corrected/";
    std::string coutput = "/Users/ka/pc/prod/hd_untiled_corrected/";
    for (Index la = from.latitude; la <= to.latitude; ++la) {
        for (Index lo = from.longitude; lo <= to.longitude; ++lo) {
            auto d__startTime = clock();
            PlanarIndex at{la, lo};
            deTile(input, output, at);
            deTile(cinput, coutput, at);
            std::cout << at << " time: " << ((double)(clock() - d__startTime) / CLOCKS_PER_SEC) << "s\n";
        }
    }
}

void deTileWithAlternative(const std::string &inputDirectory,
    const std::string &alternativeDirectory,
    const std::string &outputDirectory,
    const PlanarIndex &at) {
    if (!deTile(inputDirectory, outputDirectory, at)) {
        auto success = deTile(alternativeDirectory, outputDirectory, at);
        if (!success) {
            std::cout << "failed to detile " << at << std::endl;
        }
        return;
    }
}

void compareDEMs(const PlanarIndex &from, const PlanarIndex &to, const std::string &firstDEMDir, const std::string &secondDEMDir, const std::string &sdDEMDir) {
    TiledDEMGrid tiledGrid{firstDEMDir, from, to, {kBigTiledFileSizeX, kBigTiledFileSizeY}, {3600, 3600}, true};
    RegularDEMGrid untileddGrid{secondDEMDir, from, to, {3600, 3600}, {3600, 3600}, true};
    RegularDEMGrid sdGrid{sdDEMDir, from, to, {1200, 1200}, {1200, 1200}, true};
    std::cout << "comparing from " << from << " to " << to << endl;
    Index kDiffCutoff = 300;
    const uint8_t maxDiffs = 10;
    std::vector<PlanarIndex> diffPoints;
    Altitude maxDiff = 0;
    PlanarIndex maxDiffPlanar{0, 0};


    for (Index la = tiledGrid.indexFrom().latitude; la < tiledGrid.indexTo().latitude; ++la) {
        for (Index lo = tiledGrid.indexFrom().longitude; lo < tiledGrid.indexTo().longitude; ++lo) {
            PlanarIndex planar{la, lo};
            Altitude altDiff = std::abs(tiledGrid.get(planar) - untileddGrid.get(planar));
            if (altDiff > maxDiff) {
                maxDiff = altDiff;
                maxDiffPlanar = planar;
            }
            bool shouldContinue = false;
            for (const auto &dp : diffPoints) {
                auto diff = planar - dp;
                if (std::abs(diff.x) < 50 && std::abs(diff.y) < 50) {
                    shouldContinue = true;
                    break;;
                }
            }
            if (shouldContinue) {
                continue;
            }
            if (altDiff > kDiffCutoff) {
                diffPoints.push_back(planar);
            }
            if (diffPoints.size() == maxDiffs) {
                break;
            }
        }
        if (diffPoints.size() == maxDiffs) {
            break;
        }
    }
    std::cout << "max diff = " << maxDiff << " at " << maxDiffPlanar << " which is " << tiledGrid.planarToGeodesic(maxDiffPlanar) << endl;
    for (const auto &dp :diffPoints) {
        Altitude fa = tiledGrid.get(dp);
        Altitude sa = untileddGrid.get(dp);
        std::cout << "dp: " << dp << " " << tiledGrid.planarToGeodesic(dp) << " diff: " << fa << " - " << sa << " = " << fa - sa << endl;
    }
}

Altitude getAlt(const std::string &demsDir, Geodesic geo, bool tiled = true) {
    PlanarIndex at{Index(std::floor(geo.x)), Index(std::floor(geo.y))};
    if (tiled) {
        TiledDEMGrid region{demsDir, at, at, {kBigTiledFileSizeX, kBigTiledFileSizeY}, 3600};

        PlanarIndex planar = region.geodesicToPlanar(geo);
        return region.get(planar);
    } else {
        RegularDEMGrid region{demsDir, at, at, {3600, 3600}, 3600};
        PlanarIndex planar = region.geodesicToPlanar(geo);
        return region.get(planar);
    }
}

void applyDEMProcessing(int argc, char *argv[]) {
    if (argc > 1) {
        std::string mode(argv[1]);
        if (mode == "find_nans") {
            std::string directory = std::string(argv[2]) + "/";
            std::string output = std::string(argv[3]);
            int resolution = atoi(argv[4]);
            int diffCutoff = atoi(argv[5]);
            std::ofstream out{output};
            findNans(directory, out, resolution, diffCutoff);
        } else if (mode == "de_edem_merge_tile") {
            std::string input = std::string(argv[2]) + "/";
            std::string output = std::string(argv[3]) + "/";
            PlanarIndex at{atoi(argv[4]), atoi(argv[5])};
            deEdemMergeTile(input, output, at);
        } else if (mode == "uncompress") {
            std::string input = std::string(argv[2]);
            std::string output = std::string(argv[3]);
            uncompress(input, output);
        } else if (mode == "de_tile") {
            std::string input = std::string(argv[2]) + "/";
            std::string output = std::string(argv[3]) + "/";
            PlanarIndex at{atoi(argv[4]), atoi(argv[5])};
            deTile(input, output, at);
        } else if (mode == "de_tile_with_alternative") {
            std::string inputDir = std::string(argv[2]) + "/";
            std::string alternativeDir = std::string(argv[3]) + "/";
            std::string outputDir = std::string(argv[4]) + "/";
            PlanarIndex at{atoi(argv[5]), atoi(argv[6])};
            deTileWithAlternative(inputDir, alternativeDir, outputDir, at);
        } else if (mode == "tile") {
            std::string inputDir = std::string(argv[2]) + "/";
            std::string outputDir = std::string(argv[3]) + "/";
            PlanarIndex at{atoi(argv[4]), atoi(argv[5])};
            tile(inputDir, outputDir, at);
        } else if (mode == "sep") {
            std::string inputDir = std::string(argv[2]) + "/";
            std::string outputDir = std::string(argv[3]) + "/";
            PlanarIndex from{atoi(argv[4]), atoi(argv[5])};
            PlanarIndex to{atoi(argv[6]), atoi(argv[7])};
            prepareForCompressionNDEMFiles(inputDir, outputDir, from, to);
        } else if (mode == "slice_and_compress") {
            std::string inputDir = std::string(argv[2]) + "/";
            std::string outputDir = std::string(argv[3]) + "/";
            PlanarIndex at{atoi(argv[4]), atoi(argv[5])};
            sliceAndCompress(inputDir, outputDir, at);
        } else if (mode == "get_alt") {
            std::cout << getAlt(argv[2], {std::stod(argv[3]), std::stod(argv[4])}) << endl;
        } else {
            std::cout << "wrong mode\n";
            assert(false);
        }
    } else {
        std::cout << "no mode specified\n";
    }
}


#endif // DEM_PROCESSING_H
