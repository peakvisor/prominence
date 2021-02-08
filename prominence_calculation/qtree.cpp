#include "qtree.h"
#include "logging.hpp"

#include <algorithm>

QTree::QTree(PlanarIndex from,
             PlanarIndex to,
             PeaksContainer &peaksContainer,
             Converter converter,
             FrontierPusher frontierPusher)
             : peaksContainer(peaksContainer),
             converter(converter),
             frontierPusher(frontierPusher)
{
    std::pair<Index, Index> xMinMax =
        std::minmax(from.x, to.x);
    std::pair<Index, Index> yMinMax =
        std::minmax(from.y, to.y);

    assert(nodeResolution == pow(2, nodeResolutionLog));
    worldRealSize = {xMinMax.second - xMinMax.first,
                     yMinMax.second - yMinMax.first};
    Index maxSide = std::max(worldRealSize.x, worldRealSize.y);

    Index factoredSize = maxSide / leafSide;

    Index rootFactoredSide = 1;
    Index rootFactoredLogSide = 0;
    while (rootFactoredSide < factoredSize) {
        rootFactoredLogSide += nodeResolutionLog;
        rootFactoredSide <<= nodeResolutionLog;
    }
    
    worldLogSide = rootFactoredLogSide + leafLogSide;
    worldSide = 1 << worldLogSide;
    assert(Index(pow(2, worldLogSide)) == worldSide);
    worldOrigin = PlanarIndex{xMinMax.first, yMinMax.first};
    
    assert(worldOrigin.x == 0 && worldOrigin.y == 0);
    
    nodes.emplace_back(kQTreeIndexNone, false);
}

IslandIndex QTree::visitPoint(const Seed &seed) {
    if (seed.islandOwner == kNoIsland && !ignoreVisitTo.empty()) { // rarely
        if (ignoreVisitTo.count(seed.planar())) {
            ignoreVisitTo.erase(seed.planar());
            return kNoIsland;
        }
    }

    const PlanarIndex &visitingPlanar = seed.planar();
    
    qtree_debug(flog() << "visiting " << visitingPlanar << " island: " << seed.islandOwner << " a: " << seed.altitude << std::endl);
    
    PlanarIndex residual{visitingPlanar.x & (leafSide - 1), visitingPlanar.y & (leafSide - 1)};
    bool singleVisit = residual.x != 0 && residual.x != (leafSide - 1) && residual.y != 0 && residual.y != (leafSide - 1);
    
    mainVisit = Visit{0, worldOrigin};
    previousMainVisit = mainVisit;
    
    uint8_t visitedNeighbors = 0;
    bool onTheEdge = false;
    for (uint8_t neighborIndex = 0; neighborIndex < kMaxNeighborsCount; ++neighborIndex) {
        PlanarIndex neighborPlanar = visitingPlanar + kNeighborOffsets[neighborIndex];
        converter.normalize(neighborPlanar);

        if (!neighborPlanar.isInsideRectangleWithSize(worldOrigin, worldRealSize)) {
            neighbors[neighborIndex] = VisitingNeighbor{neighborPlanar, VisitingNeighbor::kSkipVisit, false};
            qtree_debug(flog() << "updating onTheEdge " << visitingPlanar << " from: " << neighborPlanar << " bit: " << neighborIndex << std::endl;)
            visitedNeighbors |= (1 << neighborIndex);
            onTheEdge = true;
            continue;
        }
        neighbors[neighborIndex] = VisitingNeighbor{neighborPlanar, VisitingNeighbor::kMainVisit, false};
    }
        
    int additionalVisitsCount = 0;
    Index logSide = worldLogSide;
    while (logSide > leafLogSide) {
        auto newLogSide = logSide - nodeResolutionLog;
        if (!singleVisit) {
            previousMainVisit = mainVisit;
        }
        updateVisit(mainVisit, logSide, visitingPlanar);
        if (!singleVisit) {
            for (auto &neighbor : neighbors) {
                if (neighbor.visitIndex == VisitingNeighbor::kSkipVisit) {
                    continue;
                }
                
                if (neighbor.visitIndex == VisitingNeighbor::kMainVisit) {
                    bool needNewVisit = !neighbor.planar.isInsideSquare(mainVisit.nodeOrigin, 1 << newLogSide);
                    if (needNewVisit) {
                        assert(additionalVisitsCount < 2);
                        additionalVisits[additionalVisitsCount] = previousMainVisit;
                        neighbor.visitIndex = additionalVisitsCount;
                        updateVisit(additionalVisits[additionalVisitsCount], logSide, neighbor.planar);
                        ++additionalVisitsCount;
                    }
                } else {
                    updateVisit(additionalVisits[neighbor.visitIndex], logSide, neighbor.planar);
                }
                
                debug_assert(neighbor.planar.isInsideSquare(neighbor.visitIndex == VisitingNeighbor::kMainVisit
                                           ? mainVisit.nodeOrigin
                                           : additionalVisits[neighbor.visitIndex].nodeOrigin,
                                                            1 << newLogSide));
            }
        }
        logSide = newLogSide;
    }
    assert(logSide == leafLogSide);
    
    auto &mainNode = nodes[mainVisit.node];
    bool shouldUpdate = true;
    if (seed.islandOwner != kNoIsland) {
        auto child = childIndex(mainVisit.nodeOrigin, visitingPlanar, leafLogSide, leafLogSide);
        const auto &boundaryInPlace = mainNode.getBoundaryForChild(child);
        if (boundaryInPlace.islandOwner != kNoIsland) {
            return peaksContainer.subductIslands(seed.islandOwner, boundaryInPlace.islandOwner, seed.spatial);
        }
        shouldUpdate = !mainNode.embedded.frontierSet.get(child);
    }
    
    qtree_debug(flog() << "visiting: " << visitingPlanar << " shouldUpdate: " << shouldUpdate << " singleVisit: " << singleVisit << std::endl;)
    
    IslandIndex neighborsIsland = kNoIsland;
        
    for (size_t neighborIndex = 0; neighborIndex < kMaxNeighborsCount; ++neighborIndex) {
        auto &neighbor = neighbors[neighborIndex];
        if (neighbor.visitIndex == VisitingNeighbor::kSkipVisit) {
            continue;
        }
        Visit &visit = neighbor.visitIndex == VisitingNeighbor::kMainVisit ? mainVisit : additionalVisits[neighbor.visitIndex];
        auto &neighborNode = nodes[visit.node];
        if (shouldUpdate) {
            qtree_debug(flog() << "updating " << neighbor.planar << " from visiting: " << visitingPlanar << " bit: " << (neighborIndex ^ 1) << std::endl;)
        }
        const auto &visitNeighborResults = neighborNode.visitNeighbor(visit.nodeOrigin, neighbor.planar, 1 << (neighborIndex ^ 1), shouldUpdate);
        IslandIndex nextIsland = visitNeighborResults.islandOwner;
        if (nextIsland != kNoIsland) {
            if (visitNeighborResults.shouldPushToFrontier) {
                flog() << "err: trying to update neighbors visited direction second time from visitingPlanar: " << visitingPlanar << " seed.islandOwner: " << seed.islandOwner << " i: " << neighborIndex << std::endl;
                assert(false);
            }
            if (shouldUpdate) {
                qtree_debug(flog() << "updating visiting: " << visitingPlanar << " from " << neighbor.planar << " bit: " << neighborIndex << std::endl;)
                if (visitedNeighbors & (1 << neighborIndex)) {
                    flog() << "err: trying to update visited direction second time visitingPlanar: " << visitingPlanar << " seed.islandOwner: " << seed.islandOwner << " i: " << neighborIndex << std::endl;
                } else {
                    visitedNeighbors |= (1 << neighborIndex);
                }
            }
            
            if (neighborsIsland == kNoIsland) {
                neighborsIsland = nextIsland;
                continue;
            }
            if (neighborsIsland == nextIsland) {
                continue;
            }
            neighborsIsland = peaksContainer.subductIslands(neighborsIsland, nextIsland, seed.spatial);
        } else if (visitNeighborResults.shouldPushToFrontier) {
            neighbor.shouldPushToFrontier = true;
        }
    }
    
    if (!shouldUpdate) {
        assert(seed.islandOwner != kNoIsland);
        assert(neighborsIsland != kNoIsland);
        return peaksContainer.subductIslands(seed.islandOwner, neighborsIsland, seed.spatial);
    }
    
    int island;
    if (neighborsIsland == kNoIsland) {
        assert(shouldUpdate);
        if (seed.islandOwner == kNoIsland) {
            flog() << "err:fatal: no island to set for a new boundary visited neighbors: " << visitedNeighbors << std::endl;
            assert(false);
        }
        island = seed.islandOwner;
    } else {
        if (seed.islandOwner != kNoIsland) {
//            ignoreVisitTo.insert(visitingPlanar);
            island = peaksContainer.subductIslands(seed.islandOwner, neighborsIsland, seed.spatial);
        } else {
            assert(shouldUpdate);
            island = neighborsIsland;
        }
    }
    
    island = peaksContainer.getIslandOwner(island);
    
    if (onTheEdge && peaksContainer.markUnreliableProminence(island)) {
        flog() << peaksContainer.getPeak(island).id << " prominence is unreliable becuase the island touched the edge at: " << converter.planarToGeodesic(visitingPlanar) << " a: " << seed.altitude() << std::endl;
    }
    
    if (!shouldUpdate) {
        return island;
    }
    
//    if (d__updated.count(visitingPlanar) != 0) {
//        std::cout << "double: " << visitingPlanar << std::endl;
//    } else {
//        d__updated.insert(visitingPlanar);
//    }
    
//    Seed resultSeed = seed;
//    resultSeed.islandOwner = island;
    
    for (uint8_t i = 0; i < kMaxNeighborsCount; ++i) {
        auto &neighbor = neighbors[i];
        if (neighbor.shouldPushToFrontier) {
            bool pushedToFrontier = frontierPusher(neighbor.planar);
            if (!pushedToFrontier) {
                qtree_debug(flog() << "updating " << visitingPlanar << " from not pushed to frontier: " << neighbor.planar << " bit: " << i << std::endl;)
                if (visitedNeighbors & (1 << i)) {
                    flog() << "err: trying to update neighbors visited direction from not pushed to frontier second time visitingPlanar: " << visitingPlanar << " seed.islandOwner: " << seed.islandOwner << " i: " << i << std::endl;
                }
                visitedNeighbors |= (1 << i);
            } else {
                Visit &visit = neighbor.visitIndex == VisitingNeighbor::kMainVisit ? mainVisit : additionalVisits[neighbor.visitIndex];
                auto child = childIndex(visit.nodeOrigin, neighbor.planar, leafLogSide, leafLogSide);
                if (!nodes[visit.node].embedded.frontierSet.get(child)) {
                    std::cout << "err: pushing to frontier something out of frontier set: " << neighbor.planar << " n: " << d__updated.size() << std::endl;
                    std::cout << "no frontier for child: " << child << " origin: " << visit.nodeOrigin << std::endl;
                    assert(false);
                }
                qtree_debug(flog() << "pushing: " << neighbor.planar << " to frontier from: " << visitingPlanar << std::endl;)
            }
        }
    }
    
    mainNode.update(mainVisit.nodeOrigin, visitingPlanar, island, visitedNeighbors, seed.altitude() != kOceanFloor);
    
    if (mainNode.childrenCount == 0) {
        removeNode(mainVisit.node);
    }
    
    for (int i = 0; i < additionalVisitsCount; ++i) {
        auto visitedNodeIndex = additionalVisits[i].node;
        if (nodes[visitedNodeIndex].childrenCount == 0) {
            removeNode(visitedNodeIndex);
        }
    }
    return island;
}

std::pair<QTreeIndex, PlanarIndex> QTree::nextNode(QTreeIndex nodeIndex,
                                                   PlanarIndex nodeOrigin,
                                                   Index nodeLogSide,
                                                   PlanarIndex searchedPoint,
                                                   bool addNodes)
{
    QTreeNode &node = nodes[nodeIndex];
    
    debug_assert(searchedPoint.isInsideSquare(nodeOrigin, 1 << nodeLogSide));
    PlanarIndex diff = searchedPoint - nodeOrigin;
    
    Index cellLogSide = nodeLogSide - nodeResolutionLog;
    PlanarIndex cellDiff = PlanarIndex{diff.x >> cellLogSide, diff.y >> cellLogSide};
    auto child = childIndex(nodeOrigin, searchedPoint, nodeLogSide, nodeResolutionLog);
    if (addNodes) {
        bool dontHaveThisChild = node.children[child] == kQTreeIndexNone;
        if (dontHaveThisChild) {
            Index childSide = 1 << (nodeLogSide - nodeResolutionLog);
            node.children[child] = getNewNode(nodeIndex, childSide == leafSide);
            node.childrenCount += 1;
        }
    }
    PlanarIndex newOrigin = nodeOrigin + PlanarIndex{cellDiff.x << cellLogSide, cellDiff.y << cellLogSide};
    
    debug_assert(searchedPoint.isInsideSquare(newOrigin, 1 << (nodeLogSide - nodeResolutionLog)));
    return {node.children[child], newOrigin};
}

void QTree::updateVisit(Visit &visit, const Index &logSide, const PlanarIndex &visitingPlanar) {
    debug_assert(visitingPlanar.isInsideSquare(visit.nodeOrigin, 1 << logSide));
    auto nextNodeWithOrigin = nextNode(visit.node, visit.nodeOrigin, logSide, visitingPlanar, true);
    visit.node = nextNodeWithOrigin.first;
    visit.nodeOrigin = nextNodeWithOrigin.second;
    debug_assert(visitingPlanar.isInsideSquare(visit.nodeOrigin, 1 << (logSide - nodeResolutionLog)));
}

QTreeIndex QTree::getNewNode(QTreeIndex parent, bool isLeaf) {
    QTreeIndex newNodeIndex;
    if (freeSpaces.empty()) {
        newNodeIndex = QTreeIndex(nodes.size());
        nodes.emplace_back(parent, isLeaf);
    } else {
        newNodeIndex = freeSpaces.back();
        freeSpaces.pop_back();
        nodes[newNodeIndex] = QTreeNode(parent, isLeaf);
    }
    
    return newNodeIndex;
}

void QTree::removeNode(QTreeIndex nodeIndex) {
    if (nodeIndex == 0) {
        return;
    }
    auto &node = nodes[nodeIndex];
    auto parentIndex = node.parent;
    auto &parent = nodes[parentIndex];
    qtree_debug(assert(parent.childrenCount > 0);)
    for (auto &child : parent.children) {
        if (child == nodeIndex) {
            child = kQTreeIndexNone;
            break;
        }
    }
    parent.childrenCount -= 1;
    if (parent.childrenCount == 0) {
        removeNode(parentIndex);
    }
    freeSpaces.push_back(nodeIndex);
}

void QTree::checkNonLeaf(QTreeIndex nodeIndex) {
    const auto &node = nodes[nodeIndex];
    int otherChildrenCount = 0;
    for (int i = 0; i < nodeSquare; ++i) {
        otherChildrenCount += (node.children[i] != kQTreeIndexNone);
    }
    
    qtree_debug(assert(otherChildrenCount == node.childrenCount);)
}
