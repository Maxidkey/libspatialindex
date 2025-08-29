/******************************************************************************
 * Example usage of FastRTree3D for 3D bin packing
 ******************************************************************************/

#include "FastRTree3D.h"
#include <iostream>
#include <vector>
#include <chrono>

using namespace FastRTree3D;

// Example visitor that collects intersecting cuboid IDs
class IntersectionCollector : public IVisitor {
public:
    std::vector<id_type> intersections;

    void visitData(id_type id, const Cuboid& cuboid) override {
        intersections.push_back(id);
    }

    void clear() {
        intersections.clear();
    }
};

// Example usage for 3D bin packing
int main() {
    RStarTree tree;

    // Insert some 3D cuboids (boxes in a bin packing scenario)
    std::vector<std::pair<Cuboid, id_type>> boxes;

    // Box 1: position (0,0,0) to (10,10,10)
    boxes.emplace_back(Cuboid({0.0, 0.0, 0.0}, {10.0, 10.0, 10.0}), 1);

    // Box 2: position (5,5,5) to (15,15,15) - overlaps with Box 1
    boxes.emplace_back(Cuboid({5.0, 5.0, 5.0}, {15.0, 15.0, 15.0}), 2);

    // Box 3: position (20,20,20) to (30,30,30) - no overlap
    boxes.emplace_back(Cuboid({20.0, 20.0, 20.0}, {30.0, 30.0, 30.0}), 3);

    // Box 4: position (12,12,12) to (18,18,18) - overlaps with Box 2
    boxes.emplace_back(Cuboid({12.0, 12.0, 12.0}, {18.0, 18.0, 18.0}), 4);

    // Insert boxes into R*-tree
    auto start = std::chrono::high_resolution_clock::now();

    for (const auto& box : boxes) {
        tree.insert(box.first, box.second);
    }

    auto insert_end = std::chrono::high_resolution_clock::now();

    // Query for intersections with a test cuboid
    Cuboid queryBox({6.0, 6.0, 6.0}, {14.0, 14.0, 14.0});
    IntersectionCollector collector;

    tree.intersectionQuery(queryBox, collector);

    auto query_end = std::chrono::high_resolution_clock::now();

    // Print results
    std::cout << "Inserted " << boxes.size() << " boxes into R*-tree\n";
    std::cout << "Tree height: " << tree.getHeight() << "\n";

    std::cout << "\nQuery box intersects with boxes: ";
    for (id_type id : collector.intersections) {
        std::cout << id << " ";
    }
    std::cout << "\n";

    // Performance metrics
    auto insert_time = std::chrono::duration_cast<std::chrono::microseconds>(
        insert_end - start).count();
    auto query_time = std::chrono::duration_cast<std::chrono::microseconds>(
        query_end - insert_end).count();

    std::cout << "\nPerformance:\n";
    std::cout << "Insert time: " << insert_time << " microseconds\n";
    std::cout << "Query time: " << query_time << " microseconds\n";

    // Example of checking if a new box can be placed without collision
    Cuboid newBox({1.0, 1.0, 1.0}, {9.0, 9.0, 9.0});
    collector.clear();
    tree.intersectionQuery(newBox, collector);

    if (collector.intersections.empty()) {
        std::cout << "\nNew box can be placed without collision\n";
    } else {
        std::cout << "\nNew box would collide with boxes: ";
        for (id_type id : collector.intersections) {
            std::cout << id << " ";
        }
        std::cout << "\n";
    }

    return 0;
}

// Compile with: g++ -std=c++17 -O3 FastRTree3D.cpp FastRTree3D_example.cpp -o fastrtree_example
