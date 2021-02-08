#pragma once

#include "index.h"
#include "common.h"
#include "bitmask.h"
#include <array>
#include <deque>
#include <memory>
#include <cmath>

using QTreeIndex = uint32_t;
constexpr QTreeIndex kQTreeIndexNone = -1;
constexpr int maxEmbeddedBoundary = 46;
constexpr Index nodeResolutionLog = 3;
constexpr int nodeResolution = 8;
constexpr int nodeSquare = nodeResolution * nodeResolution;
constexpr int leafLogSide = 4;
constexpr Index leafSide = 16;
constexpr Index leafSquare = leafSide * leafSide;

int childIndex(const PlanarIndex &parentOrigin, const PlanarIndex &pointInsideChild, Index parentSideLog, int parentResolutionLog);
PlanarIndex childOrigin(const PlanarIndex &parentOrigin, const Index &childSideLog, int child, int parentResolutionLog);

struct BoundaryPoint {
    static constexpr int kAllNeighborsVisited = 15;
    
    mutable IslandIndex islandOwner : 28;
    mutable unsigned int visitedNeighbors : 4;
};

struct QTreeNode {
    struct VisitNeighborResults {
        IslandIndex islandOwner;
        bool shouldPushToFrontier;
    };
    
    QTreeNode(QTreeIndex parent, bool isLeaf);
    
    void insertBoundary(int child, int islandOwner, unsigned int visitedNeighbors);
    void removeBoundary(int child, int childBoundaryIndex);
    
    VisitNeighborResults visitNeighbor(const PlanarIndex &origin, const PlanarIndex &neighbor, int neighborBit, bool shouldUpdate);
    
    void update(const PlanarIndex &origin,
                const PlanarIndex &point,
                int islandOwner,
                unsigned int visitedNeighbors,
                bool addBoundary = true);
    
    BoundaryPoint getBoundaryForChild(int child) const;

    inline bool haveAdditionalBoundary() const {
        return additionalBoundary != nullptr;
    }
    
    union {
        std::array<QTreeIndex, nodeSquare> children;
        struct {
            Bitmask<4> frontierSet;
            Bitmask<4> boundarySet;
            std::array<BoundaryPoint, maxEmbeddedBoundary> boundary;
            int boundarySize;
        } embedded;
    };
        
    QTreeIndex parent = kQTreeIndexNone;
    std::unique_ptr<std::deque<BoundaryPoint>> additionalBoundary = nullptr;
    int childrenCount = 0;
    bool isLeaf;
};
