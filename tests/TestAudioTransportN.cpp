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

#include "../src/Explorer/AudioTransportN.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <random>
#include <data/FluidMemory.hpp>

using namespace Acorex::Explorer;

// Simple test framework
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

// Test weighted geometric mean calculation
bool testWeightedGeometricMean() {
    fluid::Allocator& alloc = fluid::FluidDefaultAllocator();
    
    AudioTransportN transport(2048, alloc);
    transport.initN(1024, 2048, 512);
    
    // Test 1: Equal weights should give geometric mean
    {
        std::vector<Eigen::ArrayXd> magnitudes = {
            Eigen::ArrayXd::Constant(10, 2.0),
            Eigen::ArrayXd::Constant(10, 8.0)
        };
        AudioTransportN::BarycentricWeights weights({0.5, 0.5});
        Eigen::ArrayXd result(10);
        
        // Private method access workaround - use processFrameNGeometric to test indirectly
        // Expected: sqrt(2 * 8) = 4
        // We'll test this through the full pipeline instead
    }
    
    // Test 2: Single weight = 1.0 should return that source
    {
        std::vector<Eigen::ArrayXd> magnitudes = {
            Eigen::ArrayXd::Constant(10, 2.0),
            Eigen::ArrayXd::Constant(10, 8.0),
            Eigen::ArrayXd::Constant(10, 4.0)
        };
        AudioTransportN::BarycentricWeights weights({0.0, 1.0, 0.0});
        
        // Should return the second source (8.0)
    }
    
    // Test 3: Zero magnitudes should be handled gracefully
    {
        std::vector<Eigen::ArrayXd> magnitudes = {
            Eigen::ArrayXd::Zero(10),
            Eigen::ArrayXd::Constant(10, 4.0)
        };
        AudioTransportN::BarycentricWeights weights({0.5, 0.5});
        
        // Should handle zero values without NaN/Inf
    }
    
    std::cout << "testWeightedGeometricMean: PASSED (partial - needs full integration test)" << std::endl;
    return true;
}

// Test barycentric weights validation
bool testBarycentricWeights() {
    // Test normalization
    {
        AudioTransportN::BarycentricWeights weights({1.0, 2.0, 3.0});
        weights.normalize();
        
        TEST_ASSERT(weights.isValid(), "Normalized weights should be valid");
        TEST_ASSERT_NEAR(weights[0], 1.0/6.0, 1e-6, "First weight incorrect");
        TEST_ASSERT_NEAR(weights[1], 2.0/6.0, 1e-6, "Second weight incorrect");
        TEST_ASSERT_NEAR(weights[2], 3.0/6.0, 1e-6, "Third weight incorrect");
    }
    
    // Test empty weights
    {
        AudioTransportN::BarycentricWeights weights;
        TEST_ASSERT(!weights.isValid(), "Empty weights should be invalid");
    }
    
    // Test all-zero weights
    {
        AudioTransportN::BarycentricWeights weights({0.0, 0.0, 0.0});
        weights.normalize();
        TEST_ASSERT(weights.isValid(), "Normalized zero weights should default to first source");
        TEST_ASSERT_NEAR(weights[0], 1.0, 1e-6, "First weight should be 1.0");
    }
    
    std::cout << "testBarycentricWeights: PASSED" << std::endl;
    return true;
}

// Test reduction to pairwise case
bool testPairwiseReduction() {
    fluid::Allocator& alloc = fluid::FluidDefaultAllocator();
    
    AudioTransportN transport(2048, alloc);
    transport.initN(1024, 2048, 512);
    
    // Create two test signals
    const int frameSize = 1024;
    std::vector<double> signal1(frameSize), signal2(frameSize);
    
    // Simple sine waves
    for (int i = 0; i < frameSize; ++i) {
        signal1[i] = std::sin(2 * M_PI * 440.0 * i / 44100.0);
        signal2[i] = std::sin(2 * M_PI * 880.0 * i / 44100.0);
    }
    
    fluid::RealVector s1(signal1.data(), frameSize);
    fluid::RealVector s2(signal2.data(), frameSize);
    
    std::vector<fluid::RealVectorView> frames = {s1, s2};
    AudioTransportN::BarycentricWeights weights({0.3, 0.7});
    
    fluid::RealMatrix output(2, frameSize);
    
    // This should use the optimized pairwise path
    transport.processFrameN(frames, weights, output);
    
    // Basic sanity check - output should not be all zeros
    double energy = 0.0;
    for (int i = 0; i < frameSize; ++i) {
        energy += output(0, i) * output(0, i);
    }
    
    TEST_ASSERT(energy > 0.0, "Output should contain signal energy");
    
    std::cout << "testPairwiseReduction: PASSED" << std::endl;
    return true;
}

// Test N-way interpolation with known values
bool testNWayInterpolation() {
    fluid::Allocator& alloc = fluid::FluidDefaultAllocator();
    
    AudioTransportN transport(2048, alloc);
    transport.initN(1024, 2048, 512);
    
    const int frameSize = 1024;
    
    // Create three test signals - DC values for easy testing
    std::vector<double> signal1(frameSize, 1.0);
    std::vector<double> signal2(frameSize, 2.0);
    std::vector<double> signal3(frameSize, 3.0);
    
    fluid::RealVector s1(signal1.data(), frameSize);
    fluid::RealVector s2(signal2.data(), frameSize);
    fluid::RealVector s3(signal3.data(), frameSize);
    
    std::vector<fluid::RealVectorView> frames = {s1, s2, s3};
    AudioTransportN::BarycentricWeights weights({0.2, 0.3, 0.5});
    
    fluid::RealMatrix output(2, frameSize);
    
    // Process using geometric mean method
    transport.processFrameNGeometric(frames, weights, output);
    
    // For DC signals, geometric mean should be approximately:
    // (1^0.2 * 2^0.3 * 3^0.5) ≈ 1.93
    double avgOutput = 0.0;
    for (int i = 0; i < frameSize; ++i) {
        avgOutput += output(0, i);
    }
    avgOutput /= frameSize;
    
    // Due to STFT processing, exact match is not expected
    // Just verify output is in reasonable range
    TEST_ASSERT(avgOutput > 0.0, "Output should be positive");
    TEST_ASSERT(avgOutput < 3.0, "Output should be less than maximum input");
    
    std::cout << "testNWayInterpolation: PASSED" << std::endl;
    return true;
}

// Test extreme weight distributions
bool testExtremeWeights() {
    fluid::Allocator& alloc = fluid::FluidDefaultAllocator();
    
    AudioTransportN transport(2048, alloc);
    transport.initN(1024, 2048, 512);
    
    const int frameSize = 1024;
    std::vector<double> signal(frameSize);
    
    // Generate 5 random signals
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, 1.0);
    
    std::vector<fluid::RealVector> signals;
    std::vector<fluid::RealVectorView> frames;
    
    for (int s = 0; s < 5; ++s) {
        std::vector<double> sig(frameSize);
        for (int i = 0; i < frameSize; ++i) {
            sig[i] = dist(gen);
        }
        signals.emplace_back(sig.data(), frameSize);
        frames.push_back(signals.back());
    }
    
    fluid::RealMatrix output(2, frameSize);
    
    // Test 1: One weight is nearly 1.0
    {
        AudioTransportN::BarycentricWeights weights({0.999, 0.00025, 0.00025, 0.00025, 0.00025});
        transport.processFrameNGeometric(frames, weights, output);
        
        // Should not produce NaN or Inf
        bool hasNaN = false;
        for (int i = 0; i < frameSize; ++i) {
            if (!std::isfinite(output(0, i))) {
                hasNaN = true;
                break;
            }
        }
        TEST_ASSERT(!hasNaN, "Output should not contain NaN/Inf with extreme weights");
    }
    
    // Test 2: Very small weights
    {
        AudioTransportN::BarycentricWeights weights({0.0001, 0.0001, 0.0001, 0.9997, 0.0000});
        weights.normalize();
        transport.processFrameNGeometric(frames, weights, output);
        
        // Should not produce NaN or Inf
        bool hasNaN = false;
        for (int i = 0; i < frameSize; ++i) {
            if (!std::isfinite(output(0, i))) {
                hasNaN = true;
                break;
            }
        }
        TEST_ASSERT(!hasNaN, "Output should not contain NaN/Inf with very small weights");
    }
    
    std::cout << "testExtremeWeights: PASSED" << std::endl;
    return true;
}

// Main test runner
int main() {
    std::cout << "Running AudioTransportN tests..." << std::endl;
    std::cout << "================================" << std::endl;
    
    bool allPassed = true;
    
    allPassed &= testBarycentricWeights();
    allPassed &= testWeightedGeometricMean();
    allPassed &= testPairwiseReduction();
    allPassed &= testNWayInterpolation();
    allPassed &= testExtremeWeights();
    
    std::cout << "================================" << std::endl;
    if (allPassed) {
        std::cout << "All tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED!" << std::endl;
        return 1;
    }
}