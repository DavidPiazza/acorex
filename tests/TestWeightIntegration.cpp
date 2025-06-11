/*
The MIT License (MIT)

Copyright (c) 2024 Elowyn Fearne
*/

#include "../src/Explorer/AudioMorphEngine.h"
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <algorithms/public/KDTree.hpp>
#include <data/FluidMemory.hpp>
#include <data/FluidDataSet.hpp>
#include <data/FluidTensor.hpp>

using namespace Acorex::Explorer;

// Test framework
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "Test failed: " << message << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        return false; \
    }

// Test weight calculation integration with KD-tree query
bool testWeightCalculationIntegration() {
    // Create AudioMorphEngine
    AudioMorphEngine engine;
    engine.Initialise(44100, 512, 0.75f);
    
    // Create a simple KD-tree with 3 points in 2D space
    // First create a DataSet with proper dimensions
    fluid::FluidDataSet<std::string, double, 1> dataset(2);  // 2D points
    
    // Add points to the dataset
    fluid::RealVector point1(2);
    point1[0] = 0.0; point1[1] = 0.0;  // Point at origin
    dataset.add("0", point1);
    
    fluid::RealVector point2(2);
    point2[0] = 1.0; point2[1] = 0.0;  // Point at (1, 0)
    dataset.add("1", point2);
    
    fluid::RealVector point3(2);
    point3[0] = 0.0; point3[1] = 1.0;  // Point at (0, 1)
    dataset.add("2", point3);
    
    // Create KDTree from the dataset
    auto kdTree = std::make_shared<fluid::algorithm::KDTree>(dataset);
    
    // Set KD-tree in engine
    engine.SetKDTree(kdTree);
    
    // Test 1: Query at origin - should get exact match
    {
        fluid::RealVector query(2);
        query[0] = 0.0;
        query[1] = 0.0;
        
        auto knnResult = engine.QueryKNearest(query, 3);
        TEST_ASSERT(knnResult.first.size() == 3, "Should get 3 neighbors");
        TEST_ASSERT(knnResult.first[0] == 0.0, "First distance should be 0 (exact match)");
        
        auto weights = engine.CalculateWeights(knnResult);
        TEST_ASSERT(weights[0] == 1.0, "Exact match should get weight 1.0");
        TEST_ASSERT(weights[1] == 0.0, "Other weights should be 0");
        TEST_ASSERT(weights[2] == 0.0, "Other weights should be 0");
    }
    
    // Test 2: Query at midpoint - should get balanced weights
    {
        fluid::RealVector query(2);
        query[0] = 0.5;
        query[1] = 0.5;
        
        auto knnResult = engine.QueryKNearest(query, 3);
        TEST_ASSERT(knnResult.first.size() == 3, "Should get 3 neighbors");
        
        // All points should be at equal distance sqrt(0.5)
        double expectedDist = std::sqrt(0.5);
        TEST_ASSERT(std::abs(knnResult.first[0] - expectedDist) < 1e-6, "Distance should be sqrt(0.5)");
        TEST_ASSERT(std::abs(knnResult.first[1] - expectedDist) < 1e-6, "Distance should be sqrt(0.5)");
        
        auto weights = engine.CalculateWeights(knnResult);
        
        // With equal distances and IDW, first two weights should be equal
        TEST_ASSERT(std::abs(weights[0] - weights[1]) < 1e-6, "Equal distances should give equal weights");
    }
    
    // Test 3: Query with different weight parameters
    {
        fluid::RealVector query(2);
        query[0] = 0.2;
        query[1] = 0.1;
        
        auto knnResult = engine.QueryKNearest(query, 2);  // Only 2 neighbors
        
        // Test with different IDW exponents
        auto weights1 = engine.CalculateWeights(knnResult, false, 1.0);  // p=1
        auto weights2 = engine.CalculateWeights(knnResult, false, 3.0);  // p=3
        
        // Higher exponent should give more weight to closer point
        double ratio1 = weights1[0] / weights1[1];
        double ratio2 = weights2[0] / weights2[1];
        TEST_ASSERT(ratio2 > ratio1, "Higher exponent should increase weight ratio");
        
        // Test Gaussian kernel
        auto weightsGauss = engine.CalculateWeights(knnResult, true, 0.5);
        double sum = 0.0;
        for (size_t i = 0; i < weightsGauss.weights.size(); ++i) {
            sum += weightsGauss[i];
        }
        TEST_ASSERT(std::abs(sum - 1.0) < 1e-6, "Gaussian weights should sum to 1.0");
    }
    
    std::cout << "testWeightCalculationIntegration: PASSED" << std::endl;
    return true;
}

// Test empty KD-tree handling
bool testEmptyKDTreeHandling() {
    AudioMorphEngine engine;
    engine.Initialise(44100, 512, 0.75f);
    
    // Query without setting KD-tree
    fluid::RealVector query(2);
    query[0] = 0.0;
    query[1] = 0.0;
    
    auto knnResult = engine.QueryKNearest(query, 3);
    TEST_ASSERT(knnResult.first.size() == 0, "Empty KD-tree should return no results");
    
    auto weights = engine.CalculateWeights(knnResult);
    TEST_ASSERT(weights.weights.size() == 0, "Empty results should give empty weights");
    
    std::cout << "testEmptyKDTreeHandling: PASSED" << std::endl;
    return true;
}

// Test dimension mismatch handling
bool testDimensionMismatch() {
    AudioMorphEngine engine;
    engine.Initialise(44100, 512, 0.75f);
    
    // Create KD-tree with 3 dimensions
    fluid::FluidDataSet<std::string, double, 1> dataset(3);  // 3D points
    
    fluid::RealVector point(3);
    point[0] = 1.0; point[1] = 2.0; point[2] = 3.0;
    dataset.add("0", point);
    
    auto kdTree = std::make_shared<fluid::algorithm::KDTree>(dataset);
    
    engine.SetKDTree(kdTree);
    
    // Query with wrong dimension (2 instead of 3)
    fluid::RealVector query(2);
    query[0] = 0.0;
    query[1] = 0.0;
    
    auto knnResult = engine.QueryKNearest(query, 1);
    TEST_ASSERT(knnResult.first.size() == 0, "Dimension mismatch should return empty result");
    
    std::cout << "testDimensionMismatch: PASSED" << std::endl;
    return true;
}

int main() {
    std::cout << "Running Weight Calculation Integration tests..." << std::endl;
    std::cout << "==============================================" << std::endl;
    
    bool allPassed = true;
    
    allPassed &= testWeightCalculationIntegration();
    allPassed &= testEmptyKDTreeHandling();
    allPassed &= testDimensionMismatch();
    
    std::cout << "==============================================" << std::endl;
    if (allPassed) {
        std::cout << "All tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED!" << std::endl;
        return 1;
    }
}