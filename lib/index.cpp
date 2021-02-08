#include "index.h"
#include <array>
#include <algorithm>

PlanarIndex PlanarIndex::outputZero = {0, 0};

bool PlanarIndex::isInsideRectangle(const Rectangle &rectangle) const {
    return rectangle.origin.x <= x && x < rectangle.end.x && rectangle.origin.y <= y && y < rectangle.end.y;
}
    
std::array<PlanarIndex, kMaxFullNeighborsCount> fullImmediateVicinityTraversal(const PlanarIndex &planar) {
    static_assert(sizeof(kFullNeighborhoodTraversal) == sizeof(kFullNeighborsOffsets), "inconsistent full neighborhood traversal shifts");
    static std::array<PlanarIndex, kMaxFullNeighborsCount> traversal;
    for (unsigned i = 0; i < kMaxFullNeighborsCount; ++i) {
        traversal[i] = planar + kFullNeighborhoodTraversal[i];
    }
    return traversal;
}

std::array<PlanarIndex, kMaxNeighborsCount> neighborsScaledOffsets(int scale) {
    std::array<PlanarIndex, kMaxNeighborsCount> result;
    for (auto i = 0u; i < kMaxNeighborsCount; ++i) {
        result[i] = kNeighborOffsets[i] * scale;
    }
    return result;
}

uint8_t offsetMask(PlanarIndexDir offsetDir) {
    return 1 << offsetDir;
}

PlanarIndexDir offsetDirFromMask(uint8_t mask) {
    assert(mask != 0 && mask < (uint8_t(1) << kMaxNeighborsCount));
    mask >>= 1;
    char res = 0;
    while (mask != 0) {
        mask >>= 1;
        ++res;
    }
    return PlanarIndexDir(res);
}

void addToMask(uint8_t &mask, uint8_t offsetIndex) {
    mask |= (1 << offsetIndex);
    assert(checkMask(mask, offsetIndex));
}

bool checkMask(uint8_t mask, uint8_t offsetIndex) {
    return mask & (1 << offsetIndex);
}

PlanarIndex nextInDir(const PlanarIndex &planar, PlanarIndexDir dir) {
    assert(dir != kPlanarIndexNoDir);
    return planar + kNeighborOffsets[dir];
}
