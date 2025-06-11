/*
The MIT License (MIT)

Copyright (c) 2024 Elowyn Fearne

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "../src/Explorer/AudioMorphEngine.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <numeric>

using namespace Acorex::Explorer;

// Test framework macros
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "Test failed: " << message << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        return false; \
    }

#define TEST_ASSERT_NEAR(actual, expected, tolerance, message) \
    if (std::abs((actual) - (expected)) > (tolerance)) { \
        std::cerr << "Test failed: " << message << " at " << __FILE__ << ":" << __LINE__ \
                  << " (expected: " << expected << ", actual: " << actual << ")" << std::endl; \
        return false; \
    }

// Test inverse distance weighting with default exponent (p=2)
bool testIDWDefaultExponent() {
    AudioMorphEngine engine;
    
    // Create mock KNN result with known distances
    AudioMorphEngine::KNNResult knn;
    knn.first = {0.1, 0.2, 0.4};  // distances
    knn.second = {nullptr, nullptr, nullptr};  // IDs (not used in weight calculation)
    
    auto weights = engine.CalculateWeights(knn, false, 2.0);
    
    // Calculate expected weights manually
    // IDW with p=2: w_i = 1/d_i^2 normalized
    // w1 = 1/0.01 = 100, w2 = 1/0.04 = 25, w3 = 1/0.16 = 6.25
    // sum = 131.25
    // normalized: w1 = 100/131.25 ≈ 0.762, w2 = 25/131.25 ≈ 0.190, w3 = 6.25/131.25 ≈ 0.048
    
    TEST_ASSERT_NEAR(weights[0], 0.762, 0.001, "First weight incorrect for IDW p=2");
    TEST_ASSERT_NEAR(weights[1], 0.190, 0.001, "Second weight incorrect for IDW p=2");
    TEST_ASSERT_NEAR(weights[2], 0.048, 0.001, "Third weight incorrect for IDW p=2");
    
    // Verify weights sum to 1.0
    double sum = std::accumulate(weights.weights.begin(), weights.weights.end(), 0.0);
    TEST_ASSERT_NEAR(sum, 1.0, 1e-6, "Weights should sum to 1.0");
    
    std::cout << "testIDWDefaultExponent: PASSED" << std::endl;
    return true;
}

// Test inverse distance weighting with custom exponent
bool testIDWCustomExponent() {
    AudioMorphEngine engine;
    
    AudioMorphEngine::KNNResult knn;
    knn.first = {0.5, 1.0, 2.0};
    knn.second = {nullptr, nullptr, nullptr};
    
    // Test with p=1
    auto weights1 = engine.CalculateWeights(knn, false, 1.0);
    
    // w1 = 1/0.5 = 2, w2 = 1/1 = 1, w3 = 1/2 = 0.5
    // sum = 3.5
    // normalized: w1 = 2/3.5 ≈ 0.571, w2 = 1/3.5 ≈ 0.286, w3 = 0.5/3.5 ≈ 0.143
    
    TEST_ASSERT_NEAR(weights1[0], 0.571, 0.001, "First weight incorrect for IDW p=1");
    TEST_ASSERT_NEAR(weights1[1], 0.286, 0.001, "Second weight incorrect for IDW p=1");
    TEST_ASSERT_NEAR(weights1[2], 0.143, 0.001, "Third weight incorrect for IDW p=1");
    
    // Test with p=3
    auto weights3 = engine.CalculateWeights(knn, false, 3.0);
    
    // w1 = 1/0.125 = 8, w2 = 1/1 = 1, w3 = 1/8 = 0.125
    // sum = 9.125
    // normalized: w1 = 8/9.125 ≈ 0.877, w2 = 1/9.125 ≈ 0.110, w3 = 0.125/9.125 ≈ 0.014
    
    TEST_ASSERT_NEAR(weights3[0], 0.877, 0.001, "First weight incorrect for IDW p=3");
    TEST_ASSERT_NEAR(weights3[1], 0.110, 0.001, "Second weight incorrect for IDW p=3");
    TEST_ASSERT_NEAR(weights3[2], 0.014, 0.001, "Third weight incorrect for IDW p=3");
    
    std::cout << "testIDWCustomExponent: PASSED" << std::endl;
    return true;
}

// Test Gaussian kernel weighting
bool testGaussianKernel() {
    AudioMorphEngine engine;
    
    AudioMorphEngine::KNNResult knn;
    knn.first = {0.1, 0.3, 0.5};
    knn.second = {nullptr, nullptr, nullptr};
    
    // Test with default sigma=0.2
    auto weights = engine.CalculateWeights(knn, true, 0.2);
    
    // Gaussian: w_i = exp(-d_i^2 / (2*sigma^2))
    // sigma=0.2, 2*sigma^2 = 0.08
    // w1 = exp(-0.01/0.08) = exp(-0.125) ≈ 0.883
    // w2 = exp(-0.09/0.08) = exp(-1.125) ≈ 0.325
    // w3 = exp(-0.25/0.08) = exp(-3.125) ≈ 0.044
    
    double w1 = std::exp(-0.01/0.08);
    double w2 = std::exp(-0.09/0.08);
    double w3 = std::exp(-0.25/0.08);
    double sum = w1 + w2 + w3;
    
    TEST_ASSERT_NEAR(weights[0], w1/sum, 0.001, "First weight incorrect for Gaussian");
    TEST_ASSERT_NEAR(weights[1], w2/sum, 0.001, "Second weight incorrect for Gaussian");
    TEST_ASSERT_NEAR(weights[2], w3/sum, 0.001, "Third weight incorrect for Gaussian");
    
    // Verify weights sum to 1.0
    double weightSum = std::accumulate(weights.weights.begin(), weights.weights.end(), 0.0);
    TEST_ASSERT_NEAR(weightSum, 1.0, 1e-6, "Gaussian weights should sum to 1.0");
    
    std::cout << "testGaussianKernel: PASSED" << std::endl;
    return true;
}

// Test edge case: zero distance (exact match)
bool testZeroDistance() {
    AudioMorphEngine engine;
    
    AudioMorphEngine::KNNResult knn;
    knn.first = {0.5, 0.0, 0.3};  // Second point has zero distance
    knn.second = {nullptr, nullptr, nullptr};
    
    auto weights = engine.CalculateWeights(knn, false, 2.0);
    
    // Should return one-hot weight at index 1
    TEST_ASSERT_NEAR(weights[0], 0.0, 1e-6, "Non-zero distance should have weight 0");
    TEST_ASSERT_NEAR(weights[1], 1.0, 1e-6, "Zero distance should have weight 1");
    TEST_ASSERT_NEAR(weights[2], 0.0, 1e-6, "Non-zero distance should have weight 0");
    
    std::cout << "testZeroDistance: PASSED" << std::endl;
    return true;
}

// Test edge case: empty KNN result
bool testEmptyResult() {
    AudioMorphEngine engine;
    
    AudioMorphEngine::KNNResult knn;
    // Empty result
    
    auto weights = engine.CalculateWeights(knn, false, 2.0);
    
    // Should return empty weights (constructor default)
    TEST_ASSERT(weights.weights.size() == 0, "Empty KNN should return empty weights");
    
    std::cout << "testEmptyResult: PASSED" << std::endl;
    return true;
}

// Test edge case: single neighbor
bool testSingleNeighbor() {
    AudioMorphEngine engine;
    
    AudioMorphEngine::KNNResult knn;
    knn.first = {0.5};
    knn.second = {nullptr};
    
    auto weights = engine.CalculateWeights(knn, false, 2.0);
    
    // Single neighbor should get weight 1.0
    TEST_ASSERT(weights.weights.size() == 1, "Should have one weight");
    TEST_ASSERT_NEAR(weights[0], 1.0, 1e-6, "Single neighbor should have weight 1.0");
    
    std::cout << "testSingleNeighbor: PASSED" << std::endl;
    return true;
}

// Test very small distances (numerical stability)
bool testSmallDistances() {
    AudioMorphEngine engine;
    
    AudioMorphEngine::KNNResult knn;
    knn.first = {1e-10, 2e-10, 3e-10};  // Very small distances
    knn.second = {nullptr, nullptr, nullptr};
    
    auto weights = engine.CalculateWeights(knn, false, 2.0);
    
    // Verify weights are valid and sum to 1.0
    for (size_t i = 0; i < weights.weights.size(); ++i) {
        TEST_ASSERT(std::isfinite(weights[i]), "Weight should be finite");
        TEST_ASSERT(weights[i] >= 0.0 && weights[i] <= 1.0, "Weight should be in [0,1]");
    }
    
    double sum = std::accumulate(weights.weights.begin(), weights.weights.end(), 0.0);
    TEST_ASSERT_NEAR(sum, 1.0, 1e-6, "Weights should sum to 1.0 even with small distances");
    
    std::cout << "testSmallDistances: PASSED" << std::endl;
    return true;
}

// Test large distances
bool testLargeDistances() {
    AudioMorphEngine engine;
    
    AudioMorphEngine::KNNResult knn;
    knn.first = {100.0, 200.0, 300.0};  // Large distances
    knn.second = {nullptr, nullptr, nullptr};
    
    auto weights = engine.CalculateWeights(knn, false, 2.0);
    
    // Verify weights are valid and sum to 1.0
    for (size_t i = 0; i < weights.weights.size(); ++i) {
        TEST_ASSERT(std::isfinite(weights[i]), "Weight should be finite");
        TEST_ASSERT(weights[i] >= 0.0 && weights[i] <= 1.0, "Weight should be in [0,1]");
    }
    
    double sum = std::accumulate(weights.weights.begin(), weights.weights.end(), 0.0);
    TEST_ASSERT_NEAR(sum, 1.0, 1e-6, "Weights should sum to 1.0 even with large distances");
    
    // Verify relative ordering (closer = higher weight)
    TEST_ASSERT(weights[0] > weights[1], "Closer point should have higher weight");
    TEST_ASSERT(weights[1] > weights[2], "Closer point should have higher weight");
    
    std::cout << "testLargeDistances: PASSED" << std::endl;
    return true;
}

// Main test runner
int main() {
    std::cout << "Running Barycentric Weight Calculation tests..." << std::endl;
    std::cout << "=============================================" << std::endl;
    
    bool allPassed = true;
    
    allPassed &= testIDWDefaultExponent();
    allPassed &= testIDWCustomExponent();
    allPassed &= testGaussianKernel();
    allPassed &= testZeroDistance();
    allPassed &= testEmptyResult();
    allPassed &= testSingleNeighbor();
    allPassed &= testSmallDistances();
    allPassed &= testLargeDistances();
    
    std::cout << "=============================================" << std::endl;
    if (allPassed) {
        std::cout << "All tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED!" << std::endl;
        return 1;
    }
}