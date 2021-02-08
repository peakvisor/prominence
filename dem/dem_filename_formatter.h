#pragma once

#include "index.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

struct DEMFilenameFormatter {
    static PlanarIndex parseLatLng(const std::string &fileName) {
        std::stringstream parser{fileName};
        char northSouthIndicator, eastWestIndicator;
        Index latitude, longitude;
        parser >> northSouthIndicator >> latitude >> eastWestIndicator >> longitude;
        latitude *= northSouthIndicator == 'N' ? 1 : -1;
        longitude *= eastWestIndicator == 'E' ? 1 : -1;
        return {latitude, longitude};
    }
    
    static std::string formatLatLng(const PlanarIndex &latLng, const std::string &suffix = {}) {
        auto latitude = latLng.latitude, longitude = latLng.longitude;
        std::stringstream formatted;
        formatted << (latitude >= 0 ? 'N' : 'S');
        formatted << std::setw(2) << std::setfill('0') << std::abs(latitude);
        formatted << (longitude >= 0 ? 'E' : 'W');
        formatted << std::setw(3) << std::setfill('0') << std::abs(longitude);
        formatted << suffix;
        return formatted.str();
    }
    
    static std::string normalizedFileName(const PlanarIndex &latLng) {
        return formatLatLng(latLng, ".ndem");
    }


    static std::string ndemName(const std::string &latLng, unsigned version) {
        if (version == 0) { return latLng + ".ndem"; }
        else { return latLng + "_sanitizedV" + std::to_string(version) + ".ndem"; }
    }

    static inline std::string lastVersionPath(const Path &directory, const std::string &latLng,
        const std::function<std::string(const std::string&, unsigned)> &getVersionedFilename) {
        for (int i = 9; i >= 0; --i) {
            auto path = directory / getVersionedFilename(latLng, i);
            if (std::filesystem::exists(path)) {
                return path;
            }
        }
        return {};
    }

    static inline std::string lastVersionNdemPath(const std::string &directory, const std::string &latLng) {
        return lastVersionPath(directory, latLng, ndemName);
    }

    static inline std::string edemName(const std::string &latLng, unsigned version) {
        return "tHDDEM_c" + latLng + "_v" + std::to_string(version) + "_fM2SEP.edem";
    }

    static inline std::string lastVersionEdemPath(const std::string &directory, const std::string &latLng) {
        return lastVersionPath(directory, latLng, edemName);
    }

};
