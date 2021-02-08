#include "peak.h"
#include "distance.h"
#include "index.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>

static char kDefaultDelimiter = '\t';

void roundSix(double &w) {
    static const int shift = std::pow(10, 6);
    w = std::round(w * shift) / shift;
}

void correctIsolation(Peak &peak, const Peak &nhn) {
    if (peak.ilpLatitude == kNoDegrees ||
        DistanceTools::distanceOnEarth(peak.geo(), nhn.geo()) < DistanceTools::distanceOnEarth(peak.geo(), peak.ilp()))
    {
        peak.ilpLatitude = nhn.latitude;
        peak.ilpLongitude = nhn.longitude;
    }
}

PeakId extraPeakId(Geodesic degrees) {
    return 2e14 + PeakId((degrees.latitude + 300) * 1e4) * 1e7 + PeakId((degrees.longitude + 300) * 1e4);
}

void setProminence(Peak &peak, const Geodesic &keycolGeo, const Altitude &keycolAltitude) {
    peak.prominence = peak.altitude - keycolAltitude;
    peak.keyColLatitude = keycolGeo.latitude;
    peak.keyColLongitude = keycolGeo.longitude;
}

void setIslandParent(Peak &peak, const Peak &parentPeak, const Geodesic &keycolGeo, const Altitude &keycolAltitude) {
    setProminence(peak, keycolGeo, keycolAltitude);
    peak.islandParent = parentPeak.id;
}

std::map<PeakId, size_t> getIdToIndexMap(const std::deque<Peak> &peaks, bool resolveDuplicateIds) {
    std::map<PeakId, size_t> res;
    for (size_t i = 0; i < peaks.size(); ++i) {
        auto ins = res.try_emplace(peaks.at(i).id, i);
        if (!ins.second && resolveDuplicateIds) {
            if (peaks.at(ins.first->second).prominence <= peaks.at(i).prominence) {
                ins.first->second = i;
            }
        }
    }
    if (res.size() != peaks.size() && !resolveDuplicateIds) {
        for (size_t i = 0; i < peaks.size(); ++i) {
            if (res[peaks[i].id] != i) {
                std::cout << peaks[i].id << " is doubled\n";
            }
        }
        throw std::invalid_argument("doubled ids in getIdToIndexMap argument");
    }
    return res;
}

void sortById(std::deque<Peak> &peaks) {
    std::sort(peaks.begin(), peaks.end(), [](const Peak &f, const Peak &s) {
        return f.id < s.id;
    });
}

void checkSortedById(const std::deque<Peak> &peaks) {
    for (size_t i = 1; i < peaks.size(); ++i) {
        assert(peaks[i].id > peaks[i - 1].id);
    }
}

void dualIteration(std::deque<Peak> &first, std::deque<Peak> &second,
        const DualIterationCallback &callback,
        bool assumeSorted)
{
    if (!assumeSorted) {
        sortById(first);
        sortById(second);
    } else {
        checkSortedById(first);
        checkSortedById(second);
    }

    size_t fi = 0;
    size_t si = 0;
    size_t fsize = first.size();
    size_t ssize = second.size();
    for (;;) {
        if (fi == fsize) {
            for (; si < ssize; ++si) {
                callback(kNoIndexSizeT, si, kDualStateSecond);
            }
            break;
        }
        if (si == ssize) {
            for (; fi < fsize; ++fi) {
                callback(fi, kNoIndexSizeT, kDualStateFirst);
            }
            break;
        }
        assert(fi < fsize && si < ssize);
        const auto &fpeak = first[fi];
        const auto &speak = second[si];

        if (fpeak.id == speak.id) {
            callback(fi, si, kDualStateBoth);
            ++fi;
            ++si;
        } else if (fpeak.id < speak.id) {
            callback(fi, kNoIndexSizeT, kDualStateFirst);
            ++fi;
        } else {
            callback(kNoIndexSizeT, si, kDualStateSecond);
            ++si;
        }
    }
}

void updateSortedPeaks(std::deque<Peak> &peaks, std::deque<Peak> &updates,
        const UpdateCallback &updateCallback)
{
    dualIteration(peaks, updates, [&](size_t fi, size_t si, DualState dualState) {
        if (dualState == kDualStateBoth) {
            updateCallback(peaks[fi], updates[si]);
        }
    });
}

struct FieldHandle {
    std::string name;
    std::function<void(Peak &p, const std::string &s)> setter;
    std::function<std::string(const Peak &p)> getter;
};

template <int DecimalPoints>
double roundUpTo(double w) {
    static const int shift = std::pow(10, DecimalPoints);
    return std::round(w * shift) / shift;
}

double roundSixDoubleFromString(const std::string &s) {
    return roundUpTo<6>(std::stod(s));
}

std::string stringRoundSixFromDouble(double w) {
    return std::to_string(roundUpTo<6>(w));
}

#define FIELD_HANDLE(FIELD, FSTR, TOSTR) FieldHandle{ \
#FIELD,\
[](Peak &p, const std::string &s) { p.FIELD = FSTR(s);}, \
[](const Peak &p){ return TOSTR(p.FIELD);}}
#define DEGREES_HANDLE(FIELD) FIELD_HANDLE(FIELD, roundSixDoubleFromString, \
p.FIELD == kNoDegrees ? ""s : stringRoundSixFromDouble)
#define ID_HANDLE(FIELD) FIELD_HANDLE(FIELD, std::stoull, p.FIELD == kNoId ? ""s : std::to_string)
#define ALTITUDE_HANDLE(FIELD) FIELD_HANDLE(FIELD, std::stoi,\
p.FIELD == kNoAltitude ? ""s : std::to_string)
#define STR_HANDLE(FIELD) FIELD_HANDLE(FIELD, , )

static const std::vector<FieldHandle> kPeakFieldHandles{
    ID_HANDLE(id),
    STR_HANDLE(name),
    DEGREES_HANDLE(latitude),
    DEGREES_HANDLE(longitude),
    ALTITUDE_HANDLE(altitude),
    ALTITUDE_HANDLE(prominence),
    ID_HANDLE(islandParent),
    DEGREES_HANDLE(keyColLatitude),
    DEGREES_HANDLE(keyColLongitude),
    DEGREES_HANDLE(ilpLatitude),
    DEGREES_HANDLE(ilpLongitude),
    ID_HANDLE(nhn),
    ID_HANDLE(passId),
    ALTITUDE_HANDLE(proportionalProminence),
    FieldHandle {
        "categories",
        [](Peak &p, const std::string &s) {
            if (s.find("extraPeak"s) != std::string::npos) { p.name = ""; }
        }
    }
};

std::string Peak::toDescriptiveString() const {
    std::string res{};
    for (const auto &handle : kPeakFieldHandles) {
        if (!handle.getter) { continue; }
        auto next = handle.getter(*this);
        if (!next.empty()) {
            if (!res.empty()) { res += kDefaultDelimiter; }
            res += handle.name + ": " + next;
        }
    }
    return res;
}

std::string Peak::toStr() const {
    constexpr int kProminenceCalculatorFields = 11;
    std::string res{};
    for (int i = 0; i < kProminenceCalculatorFields; ++i) {
        res += kPeakFieldHandles.at(i).getter(*this);
        if (i != kProminenceCalculatorFields - 1) {
            res += kDefaultDelimiter;
        }
    }
    return res;
}

std::string Peak::header() {
    constexpr int kProminenceCalculatorFields = 11;
    std::string res{};
    for (int i = 0; i < kProminenceCalculatorFields; ++i) {
        res += kPeakFieldHandles.at(i).name;
        if (i != kProminenceCalculatorFields - 1) {
            res += kDefaultDelimiter;
        }
    }
    return res;
}

std::vector<std::string> split(const std::string &line, const char delimiter = '\t') {
    std::vector<std::string> res;
    static std::istringstream scanner;
    scanner.clear();
    std::string nextComponent;
    scanner.str(line);
    while (std::getline(scanner, nextComponent, delimiter)) {
        res.push_back(nextComponent);
    }
    return res;
}

std::deque<Peak> loadPeaks(const Path &inputDatabase,
        const std::function<bool(const Peak&)> &filter)
{
    std::ifstream input{inputDatabase};
    std::string header;
    std::getline(input, header);
    auto splitHeader = split(header);
    std::vector<std::optional<int>> handlesIndexes;
    for (auto &fieldName : splitHeader) {
        bool foundHandle = false;
        for (size_t i = 0; i < kPeakFieldHandles.size(); ++i) {
            if (fieldName == kPeakFieldHandles.at(i).name) {
                handlesIndexes.emplace_back(i);
                foundHandle = true;
                break;
            }
        }
        if (!foundHandle) {
            handlesIndexes.emplace_back(std::nullopt);
        }
    }

    if (!handlesIndexes.front()) {
        throw std::invalid_argument("trying to read from peaks file with invalid header: "s
            + inputDatabase.string());
    }

    std::deque<Peak> peaks;
    std::string nextLine;
    while (std::getline(input, nextLine)) {
        peaks.emplace_back();
        auto splitLine = split(nextLine);
        for (size_t i = 0; i < splitLine.size() && i < handlesIndexes.size(); ++i) {
            if (handlesIndexes.at(i) && !splitLine.at(i).empty()) {
                kPeakFieldHandles[handlesIndexes.at(i).value()].setter(
                    peaks.back(), splitLine.at(i));
            }
        }
        if (filter && !filter(peaks.back())) {
            peaks.pop_back();
        }
    }
    return peaks;
}

std::deque<Peak> loadPeaks(const Path &inputDatabase, Geodesic from, Geodesic to) {
    return loadPeaks(inputDatabase, [&from, to](const Peak &f) {
        if (f.latitude < from.latitude || f.latitude >= to.latitude) {
            return false;
        }
        return !(f.longitude < from.longitude || f.longitude >= to.longitude);
    });
}

void savePeaks(std::deque<Peak> &peaks, const Path &outputPeaksPath) {
    std::vector<bool> presentColumns(kPeakFieldHandles.size(), false);
    for (const auto &peak : peaks) {
        bool allColumnsPresent = true;
        for (size_t i = 0; i < kPeakFieldHandles.size(); ++i) {
            if (!kPeakFieldHandles[i].getter || presentColumns[i]) { continue; }
            if (!kPeakFieldHandles[i].getter(peak).empty()) { presentColumns[i] = true; }
            else { allColumnsPresent = false; }
        }
        if (allColumnsPresent) { break; }
    }

    auto lastPresentIt =
            std::find(presentColumns.rbegin(), presentColumns.rend(), true);
    auto lastPresent =
            presentColumns.size() - 1 - (lastPresentIt - presentColumns.rbegin());

    std::ofstream output(outputPeaksPath);
    for (size_t i = 0; i <= lastPresent; ++i) {
        if (presentColumns[i]) { output << kPeakFieldHandles[i].name; }
        if (i != lastPresent) { output << kDefaultDelimiter; }
        else { output << std::endl; }
    }

    for (const auto &peak : peaks) {
        for (size_t i = 0; i <= lastPresent; ++i) {
            if (presentColumns[i]) { output << kPeakFieldHandles[i].getter(peak); }
            if (i != lastPresent) { output << kDefaultDelimiter; }
            else { output << std::endl; }
        }
    }
}

#define CHECK_PEAK_FIELD_FOR_EQUALITY(FIELD) if (lhs.FIELD != rhs.FIELD) { return #FIELD;}

template <typename LhsPeakType, typename RhsPeakType>
std::string peaksAreEqual(const LhsPeakType &lhs, const RhsPeakType &rhs) {
    if (lhs.id != rhs.id) {
        throw std::invalid_argument("diff id in peaksAreEqual");
    }

    constexpr double kDistEps = 1e-3;
    if (DistanceTools::distanceOnEarth(lhs.geo(), rhs.geo()) > kDistEps) { return "geo"; }
    CHECK_PEAK_FIELD_FOR_EQUALITY(altitude)
    CHECK_PEAK_FIELD_FOR_EQUALITY(prominence)
    CHECK_PEAK_FIELD_FOR_EQUALITY(islandParent)
    CHECK_PEAK_FIELD_FOR_EQUALITY(nhn)
    CHECK_PEAK_FIELD_FOR_EQUALITY(proportionalProminence)
    CHECK_PEAK_FIELD_FOR_EQUALITY(passId)

    if (lhs.ilpLatitude != kNoDegrees) {
        bool ilpEquality = DistanceTools::distanceOnEarth(lhs.ilp(), rhs.ilp()) < kDistEps;
        if (!ilpEquality) {
            return "ilp";
        }
    }

    if (lhs.islandParent != kNoId) {
        if (DistanceTools::distanceOnEarth(lhs.keyColGeo(), rhs.keyColGeo()) > kDistEps) {
            return "keyCol";
        }
    }

    if (lhs.isExtra() != rhs.isExtra()) {
        return "namelesness";
    }

    if (!lhs.isExtra()) {
        if (lhs.name != rhs.name) {
            return "name";
        }
    }

    return "";
}

std::vector<std::pair<size_t, size_t>> getDiffPeaks(std::deque<Peak> &lhs, std::deque<Peak> &rhs) {
    std::vector<std::pair<size_t, size_t>> diff;
    auto checkEqualityLambda = [&](size_t fi, size_t si, DualState state) {
        if (state == kDualStateBoth) {
            if (!peaksAreEqual(lhs.at(fi), rhs.at(si)).empty()) {
                diff.emplace_back(fi, si);
            }
        } else {
            diff.emplace_back(fi, si);
        }
    };
    dualIteration(lhs, rhs, checkEqualityLambda);
    return diff;
}

bool peakFilesAreTheSame(const Path &lhs, const Path &rhs, bool printDiffPeaks) {
    auto lhsPeaks = loadPeaks(lhs);
    auto rhsPeaks = loadPeaks(rhs);

    auto diffs = getDiffPeaks(lhsPeaks, rhsPeaks);
    if (printDiffPeaks) {
        if (!diffs.empty()) {
            std::cout << "diff peaks:\n";
        }
        for (const auto &diff : diffs) {
            bool hasFirst = diff.first != kNoIndexSizeT;
            bool hasSecond = diff.second != kNoIndexSizeT;

            if (hasFirst && hasSecond) {
                auto diffCode = peaksAreEqual(lhsPeaks.at(diff.first), rhsPeaks.at(diff.second));
                std::cout << OUTSTR(diffCode) << std::endl;
            }

            if (hasFirst) {
                std::cout << "lhs: " << lhsPeaks.at(diff.first).toDescriptiveString() << std::endl;
            }

            if (hasSecond) {
                std::cout << "rhs: " << rhsPeaks.at(diff.second).toDescriptiveString() << std::endl;
            }

        }
    }
    return diffs.empty();
}