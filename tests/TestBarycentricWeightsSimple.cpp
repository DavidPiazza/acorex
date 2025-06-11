/*
The MIT License (MIT)

Copyright (c) 2024 Elowyn Fearne
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

// Simple test of the weight calculation algorithms without dependencies

// Barycentric weights structure
struct BarycentricWeights {
    std::vector<double> weights;
    
    BarycentricWeights() = default;
    explicit BarycentricWeights(size_t n) : weights(n, 0.0) {}
    
    double& operator[](size_t i) { return weights[i]; }
    const double& operator[](size_t i) const { return weights[i]; }
};

// KNN result structure
struct KNNResult {
    std::vector<double> first;  // distances
    std::vector<void*> second;  // IDs (not used in weight calculation)
};

// Normalize weights to sum to 1.0
void NormalizeWeights(BarycentricWeights& w) {
    double sum = std::accumulate(w.weights.begin(), w.weights.end(), 0.0);
    if (sum <= std::numeric_limits<double>::epsilon()) {
        if (!w.weights.empty()) { 
            w.weights.assign(w.weights.size(), 0.0); 
            w.weights[0] = 1.0; 
        }
        return;
    }
    for (double& v : w.weights) v /= sum;
}

// Calculate barycentric weights from KNN distances
BarycentricWeights CalculateWeights(const KNNResult& knn, bool gaussian = false, double param = 2.0) {
    BarycentricWeights weights(knn.first.size());
    size_t n = knn.first.size();
    if (n == 0) return weights;

    // Edge case: if any distance is 0, return one-hot weight
    for (size_t i = 0; i < n; ++i) {
        if (knn.first[i] == 0.0) {
            weights.weights.assign(n, 0.0);
            weights[i] = 1.0;
            return weights;
        }
    }

    std::vector<double> raw(n, 0.0);
    if (gaussian) {
        double sigma = (param > 0.0) ? param : 0.2;
        double twoSigma2 = 2.0 * sigma * sigma;
        for (size_t i = 0; i < n; ++i) {
            raw[i] = std::exp(- (knn.first[i] * knn.first[i]) / twoSigma2);
        }
    } else {
        double p = (param > 0.0) ? param : 2.0;
        for (size_t i = 0; i < n; ++i) {
            raw[i] = 1.0 / std::pow(knn.first[i], p);
        }
    }

    weights.weights = std::move(raw);
    NormalizeWeights(weights);
    return weights;
}

// Test helpers
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "Test failed: " << message << std::endl; \
        return false; \
    }

#define TEST_ASSERT_NEAR(actual, expected, tolerance, message) \
    if (std::abs((actual) - (expected)) > (tolerance)) { \
        std::cerr << "Test failed: " << message \
                  << " (expected: " << expected << ", actual: " << actual << ")" << std::endl; \
        return false; \
    }

// Test functions
bool testIDWBasic() {
    KNNResult knn;
    knn.first = {0.1, 0.2, 0.4};
    
    auto weights = CalculateWeights(knn, false, 2.0);
    
    // w1 = 1/0.01 = 100, w2 = 1/0.04 = 25, w3 = 1/0.16 = 6.25
    // sum = 131.25
    TEST_ASSERT_NEAR(weights[0], 0.762, 0.001, "IDW p=2 first weight");
    TEST_ASSERT_NEAR(weights[1], 0.190, 0.001, "IDW p=2 second weight");
    TEST_ASSERT_NEAR(weights[2], 0.048, 0.001, "IDW p=2 third weight");
    
    std::cout << "✓ IDW basic test passed" << std::endl;
    return true;
}

bool testGaussianKernel() {
    KNNResult knn;
    knn.first = {0.1, 0.3, 0.5};
    
    auto weights = CalculateWeights(knn, true, 0.2);
    
    double w1 = std::exp(-0.01/0.08);
    double w2 = std::exp(-0.09/0.08);
    double w3 = std::exp(-0.25/0.08);
    double sum = w1 + w2 + w3;
    
    TEST_ASSERT_NEAR(weights[0], w1/sum, 0.001, "Gaussian first weight");
    TEST_ASSERT_NEAR(weights[1], w2/sum, 0.001, "Gaussian second weight");
    TEST_ASSERT_NEAR(weights[2], w3/sum, 0.001, "Gaussian third weight");
    
    std::cout << "✓ Gaussian kernel test passed" << std::endl;
    return true;
}

bool testZeroDistance() {
    KNNResult knn;
    knn.first = {0.5, 0.0, 0.3};
    
    auto weights = CalculateWeights(knn, false, 2.0);
    
    TEST_ASSERT_NEAR(weights[0], 0.0, 1e-6, "Non-zero distance weight");
    TEST_ASSERT_NEAR(weights[1], 1.0, 1e-6, "Zero distance weight");
    TEST_ASSERT_NEAR(weights[2], 0.0, 1e-6, "Non-zero distance weight");
    
    std::cout << "✓ Zero distance test passed" << std::endl;
    return true;
}

bool testWeightNormalization() {
    KNNResult knn;
    knn.first = {0.1, 0.2, 0.3, 0.4, 0.5};
    
    // Test with different parameters
    for (double p = 0.5; p <= 4.0; p += 0.5) {
        auto weights = CalculateWeights(knn, false, p);
        double sum = std::accumulate(weights.weights.begin(), weights.weights.end(), 0.0);
        TEST_ASSERT_NEAR(sum, 1.0, 1e-6, "Weights should sum to 1.0");
    }
    
    std::cout << "✓ Weight normalization test passed" << std::endl;
    return true;
}

int main() {
    std::cout << "Running Barycentric Weights Tests (Simplified)" << std::endl;
    std::cout << "==============================================" << std::endl;
    
    bool allPassed = true;
    
    allPassed &= testIDWBasic();
    allPassed &= testGaussianKernel();
    allPassed &= testZeroDistance();
    allPassed &= testWeightNormalization();
    
    std::cout << "==============================================" << std::endl;
    if (allPassed) {
        std::cout << "All tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED!" << std::endl;
        return 1;
    }
}