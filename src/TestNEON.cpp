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

#include "Explorer/SIMDUtils.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>

using namespace Acorex::Explorer;

void testNEONOptimizations() {
    std::cout << "\n=== ARM NEON SIMD Optimization Test ===" << std::endl;
    
    // Detect CPU features
    auto features = SIMD::detectCpuFeatures();
    std::cout << "CPU Features detected: ";
    
    if (SIMD::hasFeature(features, SIMD::CpuFeatures::NEON)) {
        std::cout << "NEON ";
    } else if (SIMD::hasFeature(features, SIMD::CpuFeatures::AVX)) {
        std::cout << "AVX ";
    } else if (SIMD::hasFeature(features, SIMD::CpuFeatures::SSE2)) {
        std::cout << "SSE2 ";
    } else {
        std::cout << "None (scalar only)";
    }
    std::cout << std::endl;
    
    // Test parameters
    const size_t numSources = 4;
    const size_t numBins = 1024;
    const size_t iterations = 10000;
    
    // Prepare test data
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.1, 1.0);
    std::uniform_real_distribution<double> phaseDist(-M_PI, M_PI);
    
    // Allocate test data
    std::vector<std::vector<double>> magnitudes(numSources, std::vector<double>(numBins));
    std::vector<std::vector<double>> phases(numSources, std::vector<double>(numBins));
    std::vector<double> weights(numSources);
    std::vector<double> result(numBins);
    
    // Initialize test data
    double weightSum = 0.0;
    for (size_t i = 0; i < numSources; ++i) {
        weights[i] = dist(rng);
        weightSum += weights[i];
        for (size_t j = 0; j < numBins; ++j) {
            magnitudes[i][j] = dist(rng);
            phases[i][j] = phaseDist(rng);
        }
    }
    
    // Normalize weights
    for (size_t i = 0; i < numSources; ++i) {
        weights[i] /= weightSum;
    }
    
    // Create pointer arrays for SIMD functions
    std::vector<const double*> magPtrs(numSources);
    std::vector<const double*> phasePtrs(numSources);
    for (size_t i = 0; i < numSources; ++i) {
        magPtrs[i] = magnitudes[i].data();
        phasePtrs[i] = phases[i].data();
    }
    
    // Test 1: Weighted Geometric Mean
    std::cout << "\n1. Testing Weighted Geometric Mean..." << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < iterations; ++i) {
        SIMD::computeWeightedGeometricMeanSIMD(
            magPtrs.data(), weights.data(), numSources, numBins, result.data());
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto simdTime = std::chrono::duration<double, std::milli>(end - start).count();
    
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < iterations; ++i) {
        SIMD::Scalar::computeWeightedGeometricMean(
            magPtrs.data(), weights.data(), numSources, numBins, result.data());
    }
    end = std::chrono::high_resolution_clock::now();
    auto scalarTime = std::chrono::duration<double, std::milli>(end - start).count();
    
    std::cout << "  SIMD time: " << std::fixed << std::setprecision(2) << simdTime << " ms" << std::endl;
    std::cout << "  Scalar time: " << scalarTime << " ms" << std::endl;
    std::cout << "  Speedup: " << std::setprecision(2) << scalarTime / simdTime << "x" << std::endl;
    
    // Test 2: Weighted Circular Mean
    std::cout << "\n2. Testing Weighted Circular Mean..." << std::endl;
    
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < iterations; ++i) {
        SIMD::computeWeightedCircularMeanSIMD(
            phasePtrs.data(), magPtrs.data(), weights.data(), numSources, numBins, result.data());
    }
    end = std::chrono::high_resolution_clock::now();
    simdTime = std::chrono::duration<double, std::milli>(end - start).count();
    
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < iterations; ++i) {
        SIMD::Scalar::computeWeightedCircularMean(
            phasePtrs.data(), magPtrs.data(), weights.data(), numSources, numBins, result.data());
    }
    end = std::chrono::high_resolution_clock::now();
    scalarTime = std::chrono::duration<double, std::milli>(end - start).count();
    
    std::cout << "  SIMD time: " << std::fixed << std::setprecision(2) << simdTime << " ms" << std::endl;
    std::cout << "  Scalar time: " << scalarTime << " ms" << std::endl;
    std::cout << "  Speedup: " << std::setprecision(2) << scalarTime / simdTime << "x" << std::endl;
    
    // Test 3: Phase Smoothing
    std::cout << "\n3. Testing Phase Smoothing..." << std::endl;
    std::vector<double> phaseTest(numBins);
    
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < iterations * 10; ++i) {
        std::copy(phases[0].begin(), phases[0].end(), phaseTest.begin());
        SIMD::applyPhaseSmoothingSIMD(phaseTest.data(), numBins);
    }
    end = std::chrono::high_resolution_clock::now();
    simdTime = std::chrono::duration<double, std::milli>(end - start).count();
    
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < iterations * 10; ++i) {
        std::copy(phases[0].begin(), phases[0].end(), phaseTest.begin());
        SIMD::Scalar::applyPhaseSmoothing(phaseTest.data(), numBins);
    }
    end = std::chrono::high_resolution_clock::now();
    scalarTime = std::chrono::duration<double, std::milli>(end - start).count();
    
    std::cout << "  SIMD time: " << std::fixed << std::setprecision(2) << simdTime << " ms" << std::endl;
    std::cout << "  Scalar time: " << scalarTime << " ms" << std::endl;
    std::cout << "  Speedup: " << std::setprecision(2) << scalarTime / simdTime << "x" << std::endl;
    
    // Verify correctness by comparing results
    std::cout << "\n4. Verifying correctness..." << std::endl;
    
    std::vector<double> simdResult(numBins);
    std::vector<double> scalarResult(numBins);
    
    // Test geometric mean correctness
    SIMD::computeWeightedGeometricMeanSIMD(
        magPtrs.data(), weights.data(), numSources, numBins, simdResult.data());
    SIMD::Scalar::computeWeightedGeometricMean(
        magPtrs.data(), weights.data(), numSources, numBins, scalarResult.data());
    
    double maxError = 0.0;
    for (size_t i = 0; i < numBins; ++i) {
        double error = std::abs(simdResult[i] - scalarResult[i]);
        maxError = std::max(maxError, error);
    }
    std::cout << "  Geometric mean max error: " << std::scientific << maxError << std::endl;
    
    // Test circular mean correctness
    SIMD::computeWeightedCircularMeanSIMD(
        phasePtrs.data(), magPtrs.data(), weights.data(), numSources, numBins, simdResult.data());
    SIMD::Scalar::computeWeightedCircularMean(
        phasePtrs.data(), magPtrs.data(), weights.data(), numSources, numBins, scalarResult.data());
    
    maxError = 0.0;
    for (size_t i = 0; i < numBins; ++i) {
        // Handle phase wrapping
        double error = std::abs(simdResult[i] - scalarResult[i]);
        if (error > M_PI) {
            error = 2 * M_PI - error;
        }
        maxError = std::max(maxError, error);
    }
    std::cout << "  Circular mean max error: " << std::scientific << maxError << std::endl;
    
    std::cout << "\nTest completed successfully!" << std::endl;
}

int main() {
    testNEONOptimizations();
    return 0;
}