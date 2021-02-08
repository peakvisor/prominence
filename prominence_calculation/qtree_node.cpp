#include "qtree_node.h"

int childIndex(const PlanarIndex &parentOrigin, const PlanarIndex &pointInsideChild, Index parentSideLog, int parentResolutionLog) {
    debug_assert(pointInsideChild.isInsideSquare(parentOrigin, 1 << parentSideLog));
    PlanarIndex diff = pointInsideChild - parentOrigin;
    Index cellSideLog = parentSideLog - parentResolutionLog;
    PlanarIndex childDiff = PlanarIndex{diff.x >> cellSideLog, diff.y >> cellSideLog};
    
    return (childDiff.x << parentResolutionLog) + childDiff.y;
}

PlanarIndex childOrigin(const PlanarIndex &parentOrigin, const Index &childSideLog, int child, int parentResolutionLog) {
    debug_assert(child < (1 << (2 * parentResolutionLog)));
    auto x = child >> parentResolutionLog;
    auto y = child - (x << parentResolutionLog);
    PlanarIndex res{parentOrigin.x + (x << childSideLog), parentOrigin.y + (y << childSideLog)};
    debug_assert(res.isInsideSquare(parentOrigin, 1 << (childSideLog + parentResolutionLog)));
    return res;
}

QTreeNode::QTreeNode(QTreeIndex parent, bool isLeaf) : parent(parent), isLeaf(isLeaf) {
    if (isLeaf) {
        embedded.boundarySize = 0;
        embedded.frontierSet = Bitmask<4>();
        embedded.boundarySet = Bitmask<4>();
        embedded.boundary = std::array<BoundaryPoint, maxEmbeddedBoundary>();
    } else {
        children.fill(kQTreeIndexNone);
    }
}

void QTreeNode::insertBoundary(int child, int islandOwner, unsigned int visitedNeighbors) {
    assert(isLeaf);
    debug_assert(child < leafSide * leafSide);
    debug_assert(!embedded.boundarySet.get(child));
    debug_assert(!haveAdditionalBoundary()
        || (embedded.boundarySize
            == maxEmbeddedBoundary + static_cast<int>(additionalBoundary->size())));

    if (visitedNeighbors == BoundaryPoint::kAllNeighborsVisited) {
        return;
    }

    int childBoundaryIndex = embedded.boundarySet.countBefore(child);
    if (embedded.boundarySize == maxEmbeddedBoundary) {
        assert(additionalBoundary == nullptr);
        additionalBoundary = std::make_unique<std::deque<BoundaryPoint>>();
    }
    int len = -1;
    if (childBoundaryIndex < maxEmbeddedBoundary) {
        auto childIter = embedded.boundary.begin() + childBoundaryIndex;
        if (haveAdditionalBoundary()) {
            additionalBoundary->push_front(embedded.boundary.back());
        }
        len = haveAdditionalBoundary()
              ? maxEmbeddedBoundary - childBoundaryIndex - 1
              : embedded.boundarySize - childBoundaryIndex;
        std::move_backward(childIter, childIter + len, childIter + len + 1);
        embedded.boundary[childBoundaryIndex] = BoundaryPoint{islandOwner, visitedNeighbors};
    } else {
        assert(haveAdditionalBoundary());
        int childAdditionalBoundaryIndex = childBoundaryIndex - maxEmbeddedBoundary;
        additionalBoundary->insert(additionalBoundary->begin() + childAdditionalBoundaryIndex,
                                   BoundaryPoint{islandOwner, visitedNeighbors});
    }
    
    embedded.boundarySet.set(child);
    ++embedded.boundarySize;
    ++childrenCount;
}

void QTreeNode::removeBoundary(int child, int childBoundaryIndex) {
    assert(isLeaf);
    debug_assert(embedded.boundarySet.get(child));
    assert(!haveAdditionalBoundary()
        || (embedded.boundarySize
            == maxEmbeddedBoundary + static_cast<int>(additionalBoundary->size())));

    if (childrenCount == 0) {
        return;
    }
    
    embedded.boundarySet.clear(child);
    --embedded.boundarySize;
    --childrenCount;
    
    if (childBoundaryIndex >= maxEmbeddedBoundary) {
        assert(haveAdditionalBoundary());
        int childAdditionalBoundaryIndex = childBoundaryIndex - maxEmbeddedBoundary;
        additionalBoundary->erase(additionalBoundary->begin() + childAdditionalBoundaryIndex);
        if (additionalBoundary->empty()) {
            additionalBoundary = nullptr;
        }
        return;
    }
    int len = -1;
    auto childIter = embedded.boundary.begin() + childBoundaryIndex;
    if (haveAdditionalBoundary()) {
        len = maxEmbeddedBoundary - childBoundaryIndex - 1;
        std::move(childIter + 1, childIter + len + 1, childIter);
        embedded.boundary.back() = additionalBoundary->front();
        if (additionalBoundary->size() == 1) {
            additionalBoundary = nullptr;
        } else {
            additionalBoundary->pop_front();
        }
    } else {
        len = embedded.boundarySize - childBoundaryIndex;
        std::move(childIter + 1, childIter + len + 1, childIter);
    }
}

QTreeNode::VisitNeighborResults QTreeNode::visitNeighbor(const PlanarIndex &nodeOrigin,
        const PlanarIndex &neighbor, int neighborBit, bool shouldUpdate)
{
    assert(isLeaf);
    int child = childIndex(nodeOrigin, neighbor, leafLogSide, leafLogSide);
    
    if (!embedded.boundarySet.get(child)) {
        bool shouldPushToFrontier = false;
        if (shouldUpdate && !embedded.frontierSet.get(child)) {
            shouldPushToFrontier = true;
            embedded.frontierSet.set(child);
            ++childrenCount;
        }
        return {kNoIsland, shouldPushToFrontier};
    }
    debug_assert(!embedded.frontierSet.get(child));
    
    int childBoundaryIndex = embedded.boundarySet.countBefore(child);
    BoundaryPoint *neighborBoundary;
    if (childBoundaryIndex >= maxEmbeddedBoundary) {
        assert(additionalBoundary != nullptr);
        int childAdditionalBoundaryIndex = childBoundaryIndex - maxEmbeddedBoundary;
        assert(static_cast<int>(additionalBoundary->size()) > childAdditionalBoundaryIndex);
        neighborBoundary = &(*additionalBoundary)[childAdditionalBoundaryIndex];
    } else {
        neighborBoundary = &embedded.boundary[childBoundaryIndex];
    }
    
    IslandIndex island = neighborBoundary->islandOwner;
    assert(island != kNoIsland);
    assert(!(neighborBit & neighborBoundary->visitedNeighbors));
    if (shouldUpdate) {
        neighborBoundary->visitedNeighbors |= neighborBit;
        if (neighborBoundary->visitedNeighbors == BoundaryPoint::kAllNeighborsVisited) {
            removeBoundary(child, childBoundaryIndex);
        }
    }
    
    return {island, false};
}

void QTreeNode::update(const PlanarIndex &origin, const PlanarIndex &point, IslandIndex islandOwner,
        unsigned int visitedNeighbors, bool addBoundary)
{
    assert(isLeaf);
    int child = childIndex(origin, point, leafLogSide, leafLogSide);
    if (embedded.frontierSet.get(child)) {
        embedded.frontierSet.clear(child);
        --childrenCount;
    }
    if (addBoundary) {
        insertBoundary(child, islandOwner, visitedNeighbors);
    }
}

BoundaryPoint QTreeNode::getBoundaryForChild(int child) const {
    if (!embedded.boundarySet.get(child)) {
        return {kNoIsland, 0};
    }
        
    int boundaryIndex = embedded.boundarySet.countBefore(child);
    if (boundaryIndex < maxEmbeddedBoundary) {
        return embedded.boundary[boundaryIndex];
    } else {
        assert(haveAdditionalBoundary());
        return additionalBoundary->at(boundaryIndex - maxEmbeddedBoundary);
    }
}
