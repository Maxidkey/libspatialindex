/******************************************************************************
 * Project:  libsidx - C API wrapper for FastRTree3D
 * Purpose:  C API implementation for optimized 3D R*-Tree
 * Author:   Generated for FastRTree3D integration
 ******************************************************************************
 * Copyright (c) 2023
 *
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
******************************************************************************/

#include <cmath>
#include <limits>
#include <cassert>
#include <cstring>
#include <memory>
#include <vector>
#include <spatialindex/capi/sidx_api.h>
#include <spatialindex/capi/sidx_impl.h>
#include "../rtree/FastRTree3D.h"

extern "C" {

// Visitor class for collecting intersection results
class CApiVisitor : public FastRTree3D::IVisitor {
public:
    std::vector<int64_t> results;

    void visitData(FastRTree3D::id_type id, const FastRTree3D::Cuboid& cuboid) override {
        results.push_back(id);
    }

    void clear() {
        results.clear();
    }
};

// Create a new FastRTree3D instance
SIDX_DLL FastRTree3DH FastRTree3D_Create(void) {
    try {
        FastRTree3D::RStarTree* tree = new FastRTree3D::RStarTree();
        return static_cast<FastRTree3DH>(tree);
    } catch (const std::exception& e) {
        Error_PushError(RT_Failure, e.what(), "FastRTree3D_Create");
        return nullptr;
    } catch (...) {
        Error_PushError(RT_Failure, "Unknown error creating FastRTree3D", "FastRTree3D_Create");
        return nullptr;
    }
}

// Destroy a FastRTree3D instance
SIDX_DLL void FastRTree3D_Destroy(FastRTree3DH tree) {
    if (tree == nullptr) return;

    try {
        FastRTree3D::RStarTree* rtree = static_cast<FastRTree3D::RStarTree*>(tree);
        delete rtree;
    } catch (...) {
        // Ignore exceptions during destruction
    }
}

// Insert a 3D cuboid into the tree
SIDX_DLL RTError FastRTree3D_Insert(FastRTree3DH tree,
                                    int64_t id,
                                    double* pdMin,
                                    double* pdMax) {
    if (tree == nullptr) {
        Error_PushError(RT_Failure, "Tree handle is NULL", "FastRTree3D_Insert");
        return RT_Failure;
    }

    if (pdMin == nullptr || pdMax == nullptr) {
        Error_PushError(RT_Failure, "Coordinate arrays are NULL", "FastRTree3D_Insert");
        return RT_Failure;
    }

    try {
        FastRTree3D::RStarTree* rtree = static_cast<FastRTree3D::RStarTree*>(tree);

        // Create cuboid from coordinate arrays
        std::array<double, 3> low = {pdMin[0], pdMin[1], pdMin[2]};
        std::array<double, 3> high = {pdMax[0], pdMax[1], pdMax[2]};

        // Validate coordinates
        for (int i = 0; i < 3; ++i) {
            if (low[i] > high[i]) {
                Error_PushError(RT_Failure, "Invalid cuboid: min > max", "FastRTree3D_Insert");
                return RT_Failure;
            }
        }

        FastRTree3D::Cuboid cuboid(low, high);
        rtree->insert(cuboid, id);

        return RT_None;
    } catch (const std::exception& e) {
        Error_PushError(RT_Failure, e.what(), "FastRTree3D_Insert");
        return RT_Failure;
    } catch (...) {
        Error_PushError(RT_Failure, "Unknown error during insertion", "FastRTree3D_Insert");
        return RT_Failure;
    }
}

// Query for intersecting cuboids and return their IDs
SIDX_DLL RTError FastRTree3D_IntersectionQuery(FastRTree3DH tree,
                                               double* pdMin,
                                               double* pdMax,
                                               int64_t** ids,
                                               uint64_t* nResults) {
    if (tree == nullptr) {
        Error_PushError(RT_Failure, "Tree handle is NULL", "FastRTree3D_IntersectionQuery");
        return RT_Failure;
    }

    if (pdMin == nullptr || pdMax == nullptr) {
        Error_PushError(RT_Failure, "Coordinate arrays are NULL", "FastRTree3D_IntersectionQuery");
        return RT_Failure;
    }

    if (ids == nullptr || nResults == nullptr) {
        Error_PushError(RT_Failure, "Output parameters are NULL", "FastRTree3D_IntersectionQuery");
        return RT_Failure;
    }

    try {
        FastRTree3D::RStarTree* rtree = static_cast<FastRTree3D::RStarTree*>(tree);

        // Create query cuboid
        std::array<double, 3> low = {pdMin[0], pdMin[1], pdMin[2]};
        std::array<double, 3> high = {pdMax[0], pdMax[1], pdMax[2]};

        // Validate coordinates
        for (int i = 0; i < 3; ++i) {
            if (low[i] > high[i]) {
                Error_PushError(RT_Failure, "Invalid query cuboid: min > max", "FastRTree3D_IntersectionQuery");
                return RT_Failure;
            }
        }

        FastRTree3D::Cuboid queryBox(low, high);

        // Perform query
        CApiVisitor visitor;
        rtree->intersectionQuery(queryBox, visitor);

        // Allocate result array
        *nResults = visitor.results.size();
        if (*nResults > 0) {
            *ids = static_cast<int64_t*>(SIDX_NewBuffer(*nResults * sizeof(int64_t)));
            if (*ids == nullptr) {
                Error_PushError(RT_Failure, "Failed to allocate memory for results", "FastRTree3D_IntersectionQuery");
                return RT_Failure;
            }

            // Copy results
            for (uint64_t i = 0; i < *nResults; ++i) {
                (*ids)[i] = visitor.results[i];
            }
        } else {
            *ids = nullptr;
        }

        return RT_None;
    } catch (const std::exception& e) {
        Error_PushError(RT_Failure, e.what(), "FastRTree3D_IntersectionQuery");
        return RT_Failure;
    } catch (...) {
        Error_PushError(RT_Failure, "Unknown error during query", "FastRTree3D_IntersectionQuery");
        return RT_Failure;
    }
}

// Query for intersection count only (faster when you don't need the actual IDs)
SIDX_DLL RTError FastRTree3D_IntersectionQueryCount(FastRTree3DH tree,
                                                    double* pdMin,
                                                    double* pdMax,
                                                    uint64_t* nResults) {
    if (tree == nullptr) {
        Error_PushError(RT_Failure, "Tree handle is NULL", "FastRTree3D_IntersectionQueryCount");
        return RT_Failure;
    }

    if (pdMin == nullptr || pdMax == nullptr) {
        Error_PushError(RT_Failure, "Coordinate arrays are NULL", "FastRTree3D_IntersectionQueryCount");
        return RT_Failure;
    }

    if (nResults == nullptr) {
        Error_PushError(RT_Failure, "Output parameter is NULL", "FastRTree3D_IntersectionQueryCount");
        return RT_Failure;
    }

    try {
        FastRTree3D::RStarTree* rtree = static_cast<FastRTree3D::RStarTree*>(tree);

        // Create query cuboid
        std::array<double, 3> low = {pdMin[0], pdMin[1], pdMin[2]};
        std::array<double, 3> high = {pdMax[0], pdMax[1], pdMax[2]};

        // Validate coordinates
        for (int i = 0; i < 3; ++i) {
            if (low[i] > high[i]) {
                Error_PushError(RT_Failure, "Invalid query cuboid: min > max", "FastRTree3D_IntersectionQueryCount");
                return RT_Failure;
            }
        }

        FastRTree3D::Cuboid queryBox(low, high);

        // Perform query
        CApiVisitor visitor;
        rtree->intersectionQuery(queryBox, visitor);

        *nResults = visitor.results.size();

        return RT_None;
    } catch (const std::exception& e) {
        Error_PushError(RT_Failure, e.what(), "FastRTree3D_IntersectionQueryCount");
        return RT_Failure;
    } catch (...) {
        Error_PushError(RT_Failure, "Unknown error during query", "FastRTree3D_IntersectionQueryCount");
        return RT_Failure;
    }
}

// Get the height of the tree
SIDX_DLL uint32_t FastRTree3D_GetHeight(FastRTree3DH tree) {
    if (tree == nullptr) {
        Error_PushError(RT_Failure, "Tree handle is NULL", "FastRTree3D_GetHeight");
        return 0;
    }

    try {
        FastRTree3D::RStarTree* rtree = static_cast<FastRTree3D::RStarTree*>(tree);
        return rtree->getHeight();
    } catch (const std::exception& e) {
        Error_PushError(RT_Failure, e.what(), "FastRTree3D_GetHeight");
        return 0;
    } catch (...) {
        Error_PushError(RT_Failure, "Unknown error getting tree height", "FastRTree3D_GetHeight");
        return 0;
    }
}

} // extern "C"
