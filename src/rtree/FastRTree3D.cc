/******************************************************************************
 * Fast R*-Tree Implementation for 3D Cuboids - Implementation
 ******************************************************************************/

#include "FastRTree3D.h"
#include <cassert>

namespace FastRTree3D {

// Comparison functions for R*-tree splitting
static int compareLow(const void* a, const void* b, uint32_t dim) {
    const SplitEntry* e1 = static_cast<const SplitEntry*>(a);
    const SplitEntry* e2 = static_cast<const SplitEntry*>(b);

    if (e1->cuboid.low[dim] < e2->cuboid.low[dim]) return -1;
    if (e1->cuboid.low[dim] > e2->cuboid.low[dim]) return 1;
    return 0;
}

static int compareHigh(const void* a, const void* b, uint32_t dim) {
    const SplitEntry* e1 = static_cast<const SplitEntry*>(a);
    const SplitEntry* e2 = static_cast<const SplitEntry*>(b);

    if (e1->cuboid.high[dim] < e2->cuboid.high[dim]) return -1;
    if (e1->cuboid.high[dim] > e2->cuboid.high[dim]) return 1;
    return 0;
}

// LeafNode split implementation
Node* LeafNode::split(const Cuboid& cuboid, id_type id) {
    // Add the new entry temporarily
    cuboids[childCount] = cuboid;
    ids[childCount] = id;
    uint32_t totalEntries = childCount + 1;

    // R*-tree splitting algorithm
    std::vector<SplitEntry> entries;
    entries.reserve(totalEntries);

    for (uint32_t i = 0; i < totalEntries; ++i) {
        entries.emplace_back(cuboids[i], ids[i], i, 0);
    }

    uint32_t minEntries = static_cast<uint32_t>(std::ceil(totalEntries * RStarTree::SPLIT_FACTOR));
    uint32_t maxGroup1 = totalEntries - minEntries;

    double minMargin = std::numeric_limits<double>::max();
    uint32_t bestAxis = 0;
    bool bestSortOrder = false; // false = low, true = high

    // Choose split axis
    for (uint32_t axis = 0; axis < DIMENSION; ++axis) {
        // Sort by low values
        std::sort(entries.begin(), entries.end(),
                  [axis](const SplitEntry& a, const SplitEntry& b) {
                      return a.cuboid.low[axis] < b.cuboid.low[axis];
                  });

        double marginLow = 0.0;
        for (uint32_t k = minEntries; k <= maxGroup1; ++k) {
            Cuboid mbr1, mbr2;

            // First group
            mbr1 = entries[0].cuboid;
            for (uint32_t i = 1; i < k; ++i) {
                mbr1.combine(entries[i].cuboid);
            }

            // Second group
            mbr2 = entries[k].cuboid;
            for (uint32_t i = k + 1; i < totalEntries; ++i) {
                mbr2.combine(entries[i].cuboid);
            }

            marginLow += mbr1.getMargin() + mbr2.getMargin();
        }

        // Sort by high values
        std::sort(entries.begin(), entries.end(),
                  [axis](const SplitEntry& a, const SplitEntry& b) {
                      return a.cuboid.high[axis] < b.cuboid.high[axis];
                  });

        double marginHigh = 0.0;
        for (uint32_t k = minEntries; k <= maxGroup1; ++k) {
            Cuboid mbr1, mbr2;

            // First group
            mbr1 = entries[0].cuboid;
            for (uint32_t i = 1; i < k; ++i) {
                mbr1.combine(entries[i].cuboid);
            }

            // Second group
            mbr2 = entries[k].cuboid;
            for (uint32_t i = k + 1; i < totalEntries; ++i) {
                mbr2.combine(entries[i].cuboid);
            }

            marginHigh += mbr1.getMargin() + mbr2.getMargin();
        }

        double totalMargin = marginLow + marginHigh;
        if (totalMargin < minMargin) {
            minMargin = totalMargin;
            bestAxis = axis;
            bestSortOrder = (marginHigh < marginLow);
        }
    }

    // Sort entries along the chosen axis
    if (bestSortOrder) {
        std::sort(entries.begin(), entries.end(),
                  [bestAxis](const SplitEntry& a, const SplitEntry& b) {
                      return a.cuboid.high[bestAxis] < b.cuboid.high[bestAxis];
                  });
    } else {
        std::sort(entries.begin(), entries.end(),
                  [bestAxis](const SplitEntry& a, const SplitEntry& b) {
                      return a.cuboid.low[bestAxis] < b.cuboid.low[bestAxis];
                  });
    }

    // Choose split index to minimize overlap
    double minOverlap = std::numeric_limits<double>::max();
    double minArea = std::numeric_limits<double>::max();
    uint32_t bestSplit = minEntries;

    for (uint32_t k = minEntries; k <= maxGroup1; ++k) {
        Cuboid mbr1, mbr2;

        // First group
        mbr1 = entries[0].cuboid;
        for (uint32_t i = 1; i < k; ++i) {
            mbr1.combine(entries[i].cuboid);
        }

        // Second group
        mbr2 = entries[k].cuboid;
        for (uint32_t i = k + 1; i < totalEntries; ++i) {
            mbr2.combine(entries[i].cuboid);
        }

        double overlap = mbr1.getIntersectionVolume(mbr2);
        double area = mbr1.getArea() + mbr2.getArea();

        if (overlap < minOverlap || (overlap == minOverlap && area < minArea)) {
            minOverlap = overlap;
            minArea = area;
            bestSplit = k;
        }
    }

    // Create new leaf node and redistribute entries
    LeafNode* newLeaf = new LeafNode();

    // Keep first bestSplit entries in this node
    childCount = bestSplit;
    for (uint32_t i = 0; i < bestSplit; ++i) {
        cuboids[i] = entries[i].cuboid;
        ids[i] = entries[i].id;
    }

    // Move remaining entries to new node
    newLeaf->childCount = totalEntries - bestSplit;
    for (uint32_t i = bestSplit; i < totalEntries; ++i) {
        uint32_t newIndex = i - bestSplit;
        newLeaf->cuboids[newIndex] = entries[i].cuboid;
        newLeaf->ids[newIndex] = entries[i].id;
    }

    // Recompute MBRs
    recomputeMBR();
    newLeaf->recomputeMBR();

    return newLeaf;
}

// IndexNode implementation
uint32_t IndexNode::chooseSubtree(const Cuboid& cuboid) {
    if (level == 1) {
        // If pointing to leaves, use overlap minimization
        return findLeastOverlap(cuboid);
    } else {
        // Otherwise, use area enlargement minimization
        return findLeastEnlargement(cuboid);
    }
}

uint32_t IndexNode::findLeastEnlargement(const Cuboid& cuboid) {
    double minEnlargement = std::numeric_limits<double>::max();
    double minArea = std::numeric_limits<double>::max();
    uint32_t best = 0;

    for (uint32_t i = 0; i < childCount; ++i) {
        double currentArea = childMBRs[i].getArea();
        double enlargedArea = childMBRs[i].getCombined(cuboid).getArea();
        double enlargement = enlargedArea - currentArea;

        if (enlargement < minEnlargement ||
            (enlargement == minEnlargement && currentArea < minArea)) {
            minEnlargement = enlargement;
            minArea = currentArea;
            best = i;
        }
    }

    return best;
}

uint32_t IndexNode::findLeastOverlap(const Cuboid& cuboid) {
    struct OverlapEntry {
        uint32_t index;
        double enlargement;
        double area;
        double overlapIncrease;
    };

    std::vector<OverlapEntry> candidates;
    candidates.reserve(childCount);

    // Calculate enlargement for all entries
    for (uint32_t i = 0; i < childCount; ++i) {
        OverlapEntry entry;
        entry.index = i;
        entry.area = childMBRs[i].getArea();
        entry.enlargement = childMBRs[i].getCombined(cuboid).getArea() - entry.area;
        candidates.push_back(entry);
    }

    // Sort by enlargement
    std::sort(candidates.begin(), candidates.end(),
              [](const OverlapEntry& a, const OverlapEntry& b) {
                  return a.enlargement < b.enlargement;
              });

    // Consider only the first NEAR_MINIMUM_OVERLAP_FACTOR entries
    uint32_t consideredCount = std::min(childCount, NEAR_MINIMUM_OVERLAP_FACTOR);

    double minOverlapIncrease = std::numeric_limits<double>::max();
    uint32_t best = candidates[0].index;

    for (uint32_t i = 0; i < consideredCount; ++i) {
        uint32_t idx = candidates[i].index;
        Cuboid enlarged = childMBRs[idx].getCombined(cuboid);

        double overlapIncrease = 0.0;
        for (uint32_t j = 0; j < childCount; ++j) {
            if (j != idx) {
                double oldOverlap = childMBRs[idx].getIntersectionVolume(childMBRs[j]);
                double newOverlap = enlarged.getIntersectionVolume(childMBRs[j]);
                overlapIncrease += newOverlap - oldOverlap;
            }
        }

        if (overlapIncrease < minOverlapIncrease ||
            (overlapIncrease == minOverlapIncrease &&
             candidates[i].enlargement < candidates[best].enlargement)) {
            minOverlapIncrease = overlapIncrease;
            best = idx;
        }
    }

    return best;
}

void IndexNode::insert(const Cuboid& cuboid, id_type id, std::vector<Node*>& path) {
    uint32_t chosen = chooseSubtree(cuboid);
    path.push_back(this);
    children[chosen]->insert(cuboid, id, path);

    // Update MBR
    childMBRs[chosen] = children[chosen]->mbr;
    mbr.combine(cuboid);
}

Node* IndexNode::split(const Cuboid& cuboid, id_type id) {
    // This would implement index node splitting similar to leaf splitting
    // For brevity, implementing a simplified version

    // Add new entry temporarily
    uint32_t chosen = chooseSubtree(cuboid);
    Node* newChild = children[chosen]->split(cuboid, id);

    if (newChild) {
        // Need to add new child
        childMBRs[childCount] = newChild->mbr;
        children[childCount] = newChild;
        ++childCount;

        if (childCount > MAX_CAPACITY) {
            // Split this index node using similar R*-tree algorithm
            IndexNode* newIndex = new IndexNode(level);

            // Simplified split - move half the children
            uint32_t splitPoint = childCount / 2;

            newIndex->childCount = childCount - splitPoint;
            for (uint32_t i = splitPoint; i < childCount; ++i) {
                newIndex->children[i - splitPoint] = children[i];
                newIndex->childMBRs[i - splitPoint] = childMBRs[i];
            }

            childCount = splitPoint;

            recomputeMBR();
            newIndex->recomputeMBR();

            return newIndex;
        }
    }

    recomputeMBR();
    return nullptr;
}

// RStarTree implementation
void RStarTree::insert(const Cuboid& cuboid, id_type id) {
    overflowTable.assign(treeHeight, false);
    insertImpl(cuboid, id, 0); // Insert at leaf level
}

void RStarTree::insertImpl(const Cuboid& cuboid, id_type id, uint32_t targetLevel) {
    std::vector<Node*> path;

    // Find leaf node or target level node
    Node* current = root;
    while (current->level > targetLevel) {
        path.push_back(current);
        if (current->level > 0) {
            IndexNode* indexNode = static_cast<IndexNode*>(current);
            uint32_t chosen = indexNode->chooseSubtree(cuboid);
            current = indexNode->children[chosen];
        }
    }

    // Check if we can insert without overflow
    uint32_t maxCapacity = current->isLeaf() ? LeafNode::MAX_CAPACITY : IndexNode::MAX_CAPACITY;

    if (current->childCount < maxCapacity) {
        // Simple insertion
        current->insert(cuboid, id, path);
        adjustTree(path);
    } else {
        // Overflow handling
        if (!overflowTable[current->level] && !path.empty()) {
            // Try reinsertion first (R*-tree feature)
            overflowTable[current->level] = true;

            // Choose entries to reinsert
            auto reinsertIndices = chooseReinsertEntries(current, cuboid, id);

            if (!reinsertIndices.empty()) {
                // Store entries to reinsert
                std::vector<Cuboid> reinsertCuboids;
                std::vector<id_type> reinsertIds;
                std::vector<Node*> reinsertChildren; // For index nodes

                if (current->isLeaf()) {
                    LeafNode* leaf = static_cast<LeafNode*>(current);

                    // Collect entries to reinsert (in reverse order to maintain indices)
                    for (int i = reinsertIndices.size() - 1; i >= 0; --i) {
                        uint32_t idx = reinsertIndices[i];
                        reinsertCuboids.push_back(leaf->cuboids[idx]);
                        reinsertIds.push_back(leaf->ids[idx]);

                        // Remove entry by shifting remaining entries
                        for (uint32_t j = idx; j < leaf->childCount - 1; ++j) {
                            leaf->cuboids[j] = leaf->cuboids[j + 1];
                            leaf->ids[j] = leaf->ids[j + 1];
                        }
                        --leaf->childCount;
                    }

                    // Insert the new entry
                    leaf->cuboids[leaf->childCount] = cuboid;
                    leaf->ids[leaf->childCount] = id;
                    ++leaf->childCount;

                } else {
                    IndexNode* index = static_cast<IndexNode*>(current);

                    // Collect entries to reinsert
                    for (int i = reinsertIndices.size() - 1; i >= 0; --i) {
                        uint32_t idx = reinsertIndices[i];
                        reinsertCuboids.push_back(index->childMBRs[idx]);
                        reinsertIds.push_back(-1); // Index nodes don't store data IDs
                        reinsertChildren.push_back(index->children[idx]);

                        // Remove entry by shifting
                        for (uint32_t j = idx; j < index->childCount - 1; ++j) {
                            index->childMBRs[j] = index->childMBRs[j + 1];
                            index->children[j] = index->children[j + 1];
                        }
                        --index->childCount;
                    }

                    // For index nodes, we need to handle the new entry differently
                    // Since it's a split result, we'll fall back to splitting for now
                    if (reinsertChildren.empty()) {
                        // Simple case - just a regular insert
                        uint32_t chosen = index->chooseSubtree(cuboid);
                        // This shouldn't happen in normal flow, fall back to splitting
                    }
                }

                // Recompute MBR
                current->recomputeMBR();

                // Adjust tree up to root
                adjustTree(path);

                // Reinsert the removed entries
                for (size_t i = 0; i < reinsertCuboids.size(); ++i) {
                    if (current->isLeaf()) {
                        insertImpl(reinsertCuboids[i], reinsertIds[i], targetLevel);
                    } else {
                        // For index nodes, we need to reinsert the child nodes
                        // This is more complex and would require additional handling
                        // For now, fall back to splitting
                    }
                }

                return;
            }
        }

        // Fall back to splitting
        Node* newNode = current->split(cuboid, id);
        if (newNode) {
            if (path.empty()) {
                // Root split - create new root
                IndexNode* newRoot = new IndexNode(root->level + 1);
                newRoot->children[0] = root;
                newRoot->childMBRs[0] = root->mbr;
                newRoot->children[1] = newNode;
                newRoot->childMBRs[1] = newNode->mbr;
                newRoot->childCount = 2;
                newRoot->recomputeMBR();

                root = newRoot;
                ++treeHeight;
            } else {
                adjustTree(path, newNode);
            }
        } else {
            adjustTree(path);
        }
    }
}

void RStarTree::adjustTree(std::vector<Node*>& path, Node* newNode) {
    while (!path.empty()) {
        Node* parent = path.back();
        path.pop_back();

        if (newNode) {
            // Add new node to parent
            if (parent->level > 0) {
                IndexNode* indexParent = static_cast<IndexNode*>(parent);
                indexParent->childMBRs[indexParent->childCount] = newNode->mbr;
                indexParent->children[indexParent->childCount] = newNode;
                ++indexParent->childCount;

                if (indexParent->childCount > IndexNode::MAX_CAPACITY) {
                    // Parent overflow
                    newNode = indexParent->split(Cuboid(), -1);
                } else {
                    newNode = nullptr;
                }
            }
        }

        parent->recomputeMBR();

        if (path.empty() && newNode) {
            // Root split
            IndexNode* newRoot = new IndexNode(root->level + 1);
            newRoot->children[0] = root;
            newRoot->childMBRs[0] = root->mbr;
            newRoot->children[1] = newNode;
            newRoot->childMBRs[1] = newNode->mbr;
            newRoot->childCount = 2;
            newRoot->recomputeMBR();

            root = newRoot;
            ++treeHeight;
        }
    }
}

std::vector<uint32_t> RStarTree::chooseReinsertEntries(Node* node, const Cuboid& newCuboid, id_type newId) {
    // R*-tree reinsertion: remove entries that are furthest from the center
    std::vector<uint32_t> reinsertIndices;

    auto center = node->mbr.getCenter();
    uint32_t totalEntries = node->childCount + 1; // +1 for the new entry
    uint32_t reinsertCount = static_cast<uint32_t>(
        std::ceil(totalEntries * REINSERT_FACTOR));

    // Don't reinsert if it would leave too few entries
    if (reinsertCount >= totalEntries || reinsertCount == 0) {
        return reinsertIndices; // Empty vector - fall back to splitting
    }

    // Create entries with distances from center
    struct DistanceEntry {
        uint32_t index;
        double distance;
        bool isNewEntry;

        DistanceEntry(uint32_t idx, double dist, bool isNew = false)
            : index(idx), distance(dist), isNewEntry(isNew) {}
    };

    std::vector<DistanceEntry> entries;
    entries.reserve(totalEntries);

    // Calculate distances for existing entries
    if (node->isLeaf()) {
        LeafNode* leaf = static_cast<LeafNode*>(node);
        for (uint32_t i = 0; i < leaf->childCount; ++i) {
            auto entryCenter = leaf->cuboids[i].getCenter();
            double distance = 0.0;

            // Calculate squared Euclidean distance from node center
            for (uint32_t dim = 0; dim < DIMENSION; ++dim) {
                double diff = center[dim] - entryCenter[dim];
                distance += diff * diff;
            }

            entries.emplace_back(i, distance, false);
        }
    } else {
        IndexNode* index = static_cast<IndexNode*>(node);
        for (uint32_t i = 0; i < index->childCount; ++i) {
            auto entryCenter = index->childMBRs[i].getCenter();
            double distance = 0.0;

            // Calculate squared Euclidean distance from node center
            for (uint32_t dim = 0; dim < DIMENSION; ++dim) {
                double diff = center[dim] - entryCenter[dim];
                distance += diff * diff;
            }

            entries.emplace_back(i, distance, false);
        }
    }

    // Add the new entry
    auto newEntryCenter = newCuboid.getCenter();
    double newDistance = 0.0;
    for (uint32_t dim = 0; dim < DIMENSION; ++dim) {
        double diff = center[dim] - newEntryCenter[dim];
        newDistance += diff * diff;
    }
    entries.emplace_back(node->childCount, newDistance, true);

    // Sort entries by distance (furthest first)
    std::sort(entries.begin(), entries.end(),
              [](const DistanceEntry& a, const DistanceEntry& b) {
                  return a.distance > b.distance;
              });

    // Select the furthest entries for reinsertion
    // But don't reinsert the new entry - we want to insert it normally
    uint32_t collected = 0;
    for (const auto& entry : entries) {
        if (collected >= reinsertCount) break;

        // Don't reinsert the new entry
        if (!entry.isNewEntry) {
            reinsertIndices.push_back(entry.index);
            ++collected;
        }
    }

    // Sort indices in descending order for safe removal
    std::sort(reinsertIndices.begin(), reinsertIndices.end(), std::greater<uint32_t>());

    return reinsertIndices;
}

void RStarTree::intersectionQuery(const Cuboid& query, IVisitor& visitor) {
    if (root) {
        root->intersectsQuery(query, visitor);
    }
}

} // namespace FastRTree3D
