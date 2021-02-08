#include "qtree.h"
#include "distance.h"
#include <algorithm>

struct IsolationInfo {
    QTreeIndex node;
    PlanarIndex origin;
    Index logSide;
    MinMax distance;
    
    IsolationInfo(QTreeIndex node, const PlanarIndex &origin, Index logSide, const MinMax &distance)
            : node(node), origin(origin), logSide(logSide), distance(distance) {}
};

Isolation QTree::computeIsolation(const PlanarIndex &summit) const {
    std::vector<IsolationInfo> leafIsolationInfos;
    std::vector<IsolationInfo> nonLeafIsolationInfos;
    auto summitGeodesic = converter.planarToGeodesic(summit);
    auto worldOriginGeodesic = converter.planarToGeodesic(worldOrigin);
    auto worldEndGeodesic = converter.planarToGeodesic(worldOrigin + worldRealSize);
    auto distance = DistanceTools::minMaxDistanceToQuadrant(summitGeodesic, worldOriginGeodesic, worldEndGeodesic);
    
    nonLeafIsolationInfos.emplace_back(0, worldOrigin, worldLogSide, distance);
    while (!nonLeafIsolationInfos.empty()) {
        Distance minMinDistance = nonLeafIsolationInfos.front().distance.min;
        int minMinDistanceIndex = 0;
        for (int i = 1; i < static_cast<int>(nonLeafIsolationInfos.size()); ++i) {
            if (nonLeafIsolationInfos[i].distance.min < minMinDistance) {
                minMinDistance = nonLeafIsolationInfos[i].distance.min;
                minMinDistanceIndex = i;
            }
        }
        
        std::vector<IsolationInfo> nextInfos;
        auto &dividedNodeInfo = nonLeafIsolationInfos[minMinDistanceIndex];
        const QTreeNode &dividedNode = nodes[nonLeafIsolationInfos[minMinDistanceIndex].node];

        auto childLogSide = dividedNodeInfo.logSide - nodeResolutionLog;
        for (int child = 0; child < nodeSquare; ++child) {
            auto childIndex = dividedNode.children[child];
            if (childIndex == kQTreeIndexNone) {
                continue;
            }
            auto childOrig = childOrigin(dividedNodeInfo.origin, childLogSide, child, nodeResolutionLog);
            auto childSide = 1 << childLogSide;
            auto childEnd = childOrig + PlanarIndex{childSide, childSide};
            auto dist = DistanceTools::minMaxDistanceToQuadrant(summitGeodesic, converter.planarToGeodesic(childOrig), converter.planarToGeodesic(childEnd));
            nextInfos.emplace_back(childIndex, childOrig, childLogSide, dist);
        }
        
        for (int i = 0; i < static_cast<int>(nonLeafIsolationInfos.size()); ++i) {
            if (i != minMinDistanceIndex) {
                nextInfos.push_back(nonLeafIsolationInfos[i]);
            }
        }
        
        nonLeafIsolationInfos.clear();
        for (auto& info : leafIsolationInfos) {
            nextInfos.push_back(info);
        }
        
        leafIsolationInfos.clear();
        Distance minMaxDist = nextInfos.front().distance.max;
        for (auto& info : nextInfos) {
            if (info.distance.max < minMaxDist) {
                minMaxDist = info.distance.max;
            }
        }
        
        for (auto& info : nextInfos) {
            if (info.distance.min <= minMaxDist) {
                if (nodes[info.node].isLeaf) {
                    leafIsolationInfos.push_back(info);
                } else {
                    nonLeafIsolationInfos.push_back(info);
                }
            }
        }
    }
    assert(!leafIsolationInfos.empty());
    Distance isolationDistance = 4 * kEarthRadius;
    PlanarIndex isolationPoint;
    IslandIndex isolationParent = kNoIsland;
        
    for (auto& info : leafIsolationInfos) {
        const auto &node = nodes[info.node];
        for (int child = 0; child < leafSquare; ++child) {
            if (!node.embedded.boundarySet.get(child)) {
                continue;
            }
            PlanarIndex childOrig = childOrigin(info.origin, 0, child, leafLogSide);
            debug_assert(childOrig.isInsideSquare(info.origin, leafSide));
            auto to = converter.planarToGeodesic(childOrig);
            auto distance = DistanceTools::distanceOnEarth(summitGeodesic, to);
            if (distance < isolationDistance) {
                isolationDistance = distance;
                isolationPoint = childOrig;
                isolationParent = node.getBoundaryForChild(child).islandOwner;
            }
        }
    }
    
    if (isolationParent == kNoIsland) {
        flog() << "No boundary points in isolation leafs for " << summit << std::endl;
    }
    auto distanceToTheEdge = DistanceTools::distanceToEdge(summitGeodesic, worldOriginGeodesic, worldEndGeodesic);

    return Isolation{isolationDistance, isolationPoint, isolationParent, distanceToTheEdge >= isolationDistance};
}
