#pragma once

#include "index.h"
#include "bitmask.h"
#include "qtree_node.h"
#include "seed.h"
#include "peaks_container.h"

#include <deque>
#include <array>
#include <vector>
#include <functional>
#include <set>

class QTree {
private:
    struct Visit {
        QTreeIndex node;
        PlanarIndex nodeOrigin;
    };

    struct VisitingNeighbor {
        static const int kMainVisit = -1;
        static const int kSkipVisit = -2;
        
        PlanarIndex planar;
        int visitIndex;
        bool shouldPushToFrontier;
    };
    
public:
    using FrontierPusher = std::function<bool(const PlanarIndex&)>;
    
    QTree(PlanarIndex from,
          PlanarIndex to,
          PeaksContainer &peaksContainer,
          Converter converter,
          FrontierPusher frontierPusher);
    
    Isolation computeIsolation(const PlanarIndex &summit) const;
    
    IslandIndex visitPoint(const Seed &seed);
    
    size_t nodesSize() const {
        return nodes.size();
    }
    
private:
    std::pair<QTreeIndex, PlanarIndex> nextNode(QTreeIndex nodeIndex,
                                                PlanarIndex nodeOrigin,
                                                Index nodeSide,
                                                PlanarIndex searchedPoint,
                                                bool addNodes);
    void updateVisit(Visit &visit, const Index &logSide, const PlanarIndex& visitingPlanar);
    QTreeIndex getNewNode(QTreeIndex parent, bool isLeaf);
    void removeNode(QTreeIndex nodeIndex);
    void checkNonLeaf(QTreeIndex nodeIndex);
    
    PlanarIndex worldOrigin;
    PlanarIndex worldRealSize;
    Index worldSide;
    Index worldLogSide;
    Index maxIndex;
    Index minIndex;
    
    PeaksContainer &peaksContainer;
    Converter converter;
    FrontierPusher frontierPusher;
    
    Visit mainVisit;
    Visit previousMainVisit;
    std::array<Visit, 2> additionalVisits;
    std::array<VisitingNeighbor, 4> neighbors;
    
    std::deque<QTreeNode> nodes;
    std::vector<QTreeIndex> freeSpaces;
    
    std::set<PlanarIndex> ignoreVisitTo;
    std::set<PlanarIndex> d__updated;
    
    uint64_t d__visited = 0;
};
