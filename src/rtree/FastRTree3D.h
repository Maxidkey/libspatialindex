/******************************************************************************
 * Fast R*-Tree Implementation for 3D Cuboids
 * Optimized for 3D bin packing algorithms
 *
 * Features:
 * - Fixed 3D dimension (no dynamic dimensionality)
 * - No I/O operations (memory-only)
 * - No data storage in leaves (ID-only)
 * - R*-tree splitting algorithm only
 * - Minimal virtual function calls
 * - Stack-based memory allocation where possible
 ******************************************************************************/

#pragma once

#include <vector>
#include <stack>
#include <array>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cstring>

namespace FastRTree3D {

using id_type = int64_t;
constexpr uint32_t DIMENSION = 3;

// 3D Cuboid representation
struct Cuboid {
    std::array<double, DIMENSION> low;
    std::array<double, DIMENSION> high;

    Cuboid() {
        for (int i = 0; i < DIMENSION; ++i) {
            low[i] = std::numeric_limits<double>::max();
            high[i] = -std::numeric_limits<double>::max();
        }
    }

    Cuboid(const std::array<double, DIMENSION>& l, const std::array<double, DIMENSION>& h)
        : low(l), high(h) {}

    double getArea() const {
        double area = 1.0;
        for (int i = 0; i < DIMENSION; ++i) {
            area *= (high[i] - low[i]);
        }
        return area;
    }

    double getMargin() const {
        double margin = 0.0;
        for (int i = 0; i < DIMENSION; ++i) {
            margin += (high[i] - low[i]);
        }
        return 2.0 * margin;
    }

    bool intersects(const Cuboid& other) const {
        for (int i = 0; i < DIMENSION; ++i) {
            if (low[i] > other.high[i] || high[i] < other.low[i]) {
                return false;
            }
        }
        return true;
    }

    bool contains(const Cuboid& other) const {
        for (int i = 0; i < DIMENSION; ++i) {
            if (low[i] > other.low[i] || high[i] < other.high[i]) {
                return false;
            }
        }
        return true;
    }

    void combine(const Cuboid& other) {
        for (int i = 0; i < DIMENSION; ++i) {
            low[i] = std::min(low[i], other.low[i]);
            high[i] = std::max(high[i], other.high[i]);
        }
    }

    Cuboid getCombined(const Cuboid& other) const {
        Cuboid result = *this;
        result.combine(other);
        return result;
    }

    double getIntersectionVolume(const Cuboid& other) const {
        double volume = 1.0;
        for (int i = 0; i < DIMENSION; ++i) {
            double overlap = std::min(high[i], other.high[i]) - std::max(low[i], other.low[i]);
            if (overlap <= 0.0) return 0.0;
            volume *= overlap;
        }
        return volume;
    }

    std::array<double, DIMENSION> getCenter() const {
        std::array<double, DIMENSION> center;
        for (int i = 0; i < DIMENSION; ++i) {
            center[i] = (low[i] + high[i]) * 0.5;
        }
        return center;
    }
};

// Forward declarations
class Node;
class LeafNode;
class IndexNode;

// Visitor interface for queries
class IVisitor {
public:
    virtual ~IVisitor() = default;
    virtual void visitData(id_type id, const Cuboid& cuboid) = 0;
};

// Entry for R*-tree splitting
struct SplitEntry {
    Cuboid cuboid;
    id_type id;
    uint32_t index;
    uint32_t sortDim;

    SplitEntry(const Cuboid& c, id_type i, uint32_t idx, uint32_t dim)
        : cuboid(c), id(i), index(idx), sortDim(dim) {}
};

// Base node class
class Node {
public:
    uint32_t level;
    uint32_t childCount;
    Cuboid mbr;

    Node(uint32_t lvl) : level(lvl), childCount(0) {}
    virtual ~Node() = default;

    bool isLeaf() const { return level == 0; }
    virtual bool intersectsQuery(const Cuboid& query, IVisitor& visitor) = 0;
    virtual void insert(const Cuboid& cuboid, id_type id, std::vector<Node*>& path) = 0;
    virtual Node* split(const Cuboid& cuboid, id_type id) = 0;
    virtual void recomputeMBR() = 0;
};

// Leaf node implementation
class LeafNode : public Node {
public:
    static constexpr uint32_t MAX_CAPACITY = 16; // Tunable parameter

    std::array<Cuboid, MAX_CAPACITY + 1> cuboids; // +1 for overflow during split
    std::array<id_type, MAX_CAPACITY + 1> ids;

    LeafNode() : Node(0) {}

    bool intersectsQuery(const Cuboid& query, IVisitor& visitor) override {
        bool found = false;
        for (uint32_t i = 0; i < childCount; ++i) {
            if (cuboids[i].intersects(query)) {
                visitor.visitData(ids[i], cuboids[i]);
                found = true;
            }
        }
        return found;
    }

    void insert(const Cuboid& cuboid, id_type id, std::vector<Node*>& path) override {
        cuboids[childCount] = cuboid;
        ids[childCount] = id;
        ++childCount;
        mbr.combine(cuboid);
    }

    Node* split(const Cuboid& cuboid, id_type id) override;

    void recomputeMBR() override {
        if (childCount == 0) return;

        mbr = cuboids[0];
        for (uint32_t i = 1; i < childCount; ++i) {
            mbr.combine(cuboids[i]);
        }
    }
};

// Index node implementation
class IndexNode : public Node {
public:
    static constexpr uint32_t MAX_CAPACITY = 16; // Tunable parameter

    std::array<Cuboid, MAX_CAPACITY + 1> childMBRs; // +1 for overflow during split
    std::array<Node*, MAX_CAPACITY + 1> children;

    IndexNode(uint32_t lvl) : Node(lvl) {}

    ~IndexNode() override {
        for (uint32_t i = 0; i < childCount; ++i) {
            delete children[i];
        }
    }

    bool intersectsQuery(const Cuboid& query, IVisitor& visitor) override {
        bool found = false;
        for (uint32_t i = 0; i < childCount; ++i) {
            if (childMBRs[i].intersects(query)) {
                if (children[i]->intersectsQuery(query, visitor)) {
                    found = true;
                }
            }
        }
        return found;
    }

    void insert(const Cuboid& cuboid, id_type id, std::vector<Node*>& path) override;
    Node* split(const Cuboid& cuboid, id_type id) override;

    void recomputeMBR() override {
        if (childCount == 0) return;

        mbr = childMBRs[0];
        for (uint32_t i = 1; i < childCount; ++i) {
            mbr.combine(childMBRs[i]);
        }
    }

private:
    uint32_t chooseSubtree(const Cuboid& cuboid);
    uint32_t findLeastOverlap(const Cuboid& cuboid);
    uint32_t findLeastEnlargement(const Cuboid& cuboid);
};

// Main R*-Tree class
class RStarTree {
public:
    static constexpr double REINSERT_FACTOR = 0.3;
    static constexpr double SPLIT_FACTOR = 0.4;
    static constexpr uint32_t NEAR_MINIMUM_OVERLAP_FACTOR = 32;

private:
    Node* root;
    uint32_t treeHeight;
    std::vector<bool> overflowTable; // Per-level overflow tracking

public:
    RStarTree() : root(new LeafNode()), treeHeight(1) {}

    ~RStarTree() {
        delete root;
    }

    // Non-copyable for simplicity
    RStarTree(const RStarTree&) = delete;
    RStarTree& operator=(const RStarTree&) = delete;

    void insert(const Cuboid& cuboid, id_type id);
    void intersectionQuery(const Cuboid& query, IVisitor& visitor);

    uint32_t getHeight() const { return treeHeight; }

private:
    void insertImpl(const Cuboid& cuboid, id_type id, uint32_t targetLevel);
    void adjustTree(std::vector<Node*>& path, Node* newNode = nullptr);
    std::vector<uint32_t> chooseReinsertEntries(Node* node, const Cuboid& newCuboid, id_type newId);
};

} // namespace FastRTree3D
