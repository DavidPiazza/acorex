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

#include "SIMDUtils.hpp"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <memory>

namespace Acorex {
namespace Explorer {

class SIMDBenchmark {
public:
    struct BenchmarkResult {
        std::string name;
        double scalarTime;
        double simdTime;
        double speedup;
        size_t iterations;
    };
    
    static void runAllBenchmarks() {
        std::cout << "\n=== SIMD Performance Benchmarks ===" << std::endl;
        std::cout << "CPU Features detected: " << getCpuFeaturesString() << std::endl;
        std::cout << std::endl;
        
        std::vector<BenchmarkResult> results;
        
        // Benchmark configurations
        const std::vector<size_t> binCounts = {512, 1024, 2048, 4096};
        const std::vector<size_t> sourceCounts = {2, 4, 6, 8};
        const size_t iterations = 1000;
        
        // Run benchmarks for different configurations
        for (size_t bins : binCounts) {
            for (size_t sources : sourceCounts) {
                results.push_back(benchmarkWeightedGeometricMean(sources, bins, iterations));
                results.push_back(benchmarkWeightedCircularMean(sources, bins, iterations));
            }
        }
        
        // Add phase smoothing benchmark
        for (size_t bins : binCounts) {
            results.push_back(benchmarkPhaseSmoothing(bins, iterations * 10));
        }
        
        // Print results table
        printResults(results);
    }
    
private:
    static std::string getCpuFeaturesString() {
        auto features = SIMD::detectCpuFeatures();
        std::string result;
        
        if (SIMD::hasFeature(features, SIMD::CpuFeatures::SSE2)) result += "SSE2 ";
        if (SIMD::hasFeature(features, SIMD::CpuFeatures::SSE41)) result += "SSE4.1 ";
        if (SIMD::hasFeature(features, SIMD::CpuFeatures::AVX)) result += "AVX ";
        if (SIMD::hasFeature(features, SIMD::CpuFeatures::AVX2)) result += "AVX2 ";
        if (SIMD::hasFeature(features, SIMD::CpuFeatures::FMA)) result += "FMA ";
        if (SIMD::hasFeature(features, SIMD::CpuFeatures::NEON)) result += "NEON ";
        
        if (result.empty()) result = "None";
        return result;
    }
    
    static BenchmarkResult benchmarkWeightedGeometricMean(
        size_t numSources, size_t numBins, size_t iterations) {
        
        // Prepare test data
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        
        // Allocate aligned memory for better SIMD performance
        std::vector<std::unique_ptr<double[]>> magnitudes(numSources);
        std::vector<const double*> magPtrs(numSources);
        
        for (size_t i = 0; i < numSources; ++i) {
            magnitudes[i].reset(new double[numBins]);
            magPtrs[i] = magnitudes[i].get();
            
            for (size_t j = 0; j < numBins; ++j) {
                magnitudes[i][j] = dist(rng);
            }
        }
        
        // Generate barycentric weights
        std::vector<double> weights(numSources);
        double totalWeight = 0.0;
        for (size_t i = 0; i < numSources; ++i) {
            weights[i] = dist(rng);
            totalWeight += weights[i];
        }
        for (auto& w : weights) w /= totalWeight;
        
        std::unique_ptr<double[]> result(new double[numBins]);
        
        // Benchmark scalar implementation
        auto startScalar = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < iterations; ++iter) {
            SIMD::Scalar::computeWeightedGeometricMean(
                magPtrs.data(), weights.data(), numSources, numBins, result.get());
        }
        auto endScalar = std::chrono::high_resolution_clock::now();
        
        // Benchmark SIMD implementation
        auto startSIMD = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < iterations; ++iter) {
            SIMD::computeWeightedGeometricMeanSIMD(
                magPtrs.data(), weights.data(), numSources, numBins, result.get());
        }
        auto endSIMD = std::chrono::high_resolution_clock::now();
        
        // Calculate times
        double scalarTime = std::chrono::duration<double, std::milli>(
            endScalar - startScalar).count();
        double simdTime = std::chrono::duration<double, std::milli>(
            endSIMD - startSIMD).count();
        
        return {
            "Weighted Geometric Mean (" + std::to_string(numSources) + " sources, " + 
            std::to_string(numBins) + " bins)",
            scalarTime,
            simdTime,
            scalarTime / simdTime,
            iterations
        };
    }
    
    static BenchmarkResult benchmarkWeightedCircularMean(
        size_t numSources, size_t numBins, size_t iterations) {
        
        // Prepare test data
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> distMag(0.0, 1.0);
        std::uniform_real_distribution<double> distPhase(-M_PI, M_PI);
        
        std::vector<std::unique_ptr<double[]>> magnitudes(numSources);
        std::vector<std::unique_ptr<double[]>> phases(numSources);
        std::vector<const double*> magPtrs(numSources);
        std::vector<const double*> phasePtrs(numSources);
        
        for (size_t i = 0; i < numSources; ++i) {
            magnitudes[i].reset(new double[numBins]);
            phases[i].reset(new double[numBins]);
            magPtrs[i] = magnitudes[i].get();
            phasePtrs[i] = phases[i].get();
            
            for (size_t j = 0; j < numBins; ++j) {
                magnitudes[i][j] = distMag(rng);
                phases[i][j] = distPhase(rng);
            }
        }
        
        // Generate barycentric weights
        std::vector<double> weights(numSources);
        double totalWeight = 0.0;
        for (size_t i = 0; i < numSources; ++i) {
            weights[i] = distMag(rng);
            totalWeight += weights[i];
        }
        for (auto& w : weights) w /= totalWeight;
        
        std::unique_ptr<double[]> result(new double[numBins]);
        
        // Benchmark scalar implementation
        auto startScalar = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < iterations; ++iter) {
            SIMD::Scalar::computeWeightedCircularMean(
                phasePtrs.data(), magPtrs.data(), weights.data(), 
                numSources, numBins, result.get());
        }
        auto endScalar = std::chrono::high_resolution_clock::now();
        
        // Benchmark SIMD implementation
        auto startSIMD = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < iterations; ++iter) {
            SIMD::computeWeightedCircularMeanSIMD(
                phasePtrs.data(), magPtrs.data(), weights.data(), 
                numSources, numBins, result.get());
        }
        auto endSIMD = std::chrono::high_resolution_clock::now();
        
        // Calculate times
        double scalarTime = std::chrono::duration<double, std::milli>(
            endScalar - startScalar).count();
        double simdTime = std::chrono::duration<double, std::milli>(
            endSIMD - startSIMD).count();
        
        return {
            "Weighted Circular Mean (" + std::to_string(numSources) + " sources, " + 
            std::to_string(numBins) + " bins)",
            scalarTime,
            simdTime,
            scalarTime / simdTime,
            iterations
        };
    }
    
    static BenchmarkResult benchmarkPhaseSmoothing(size_t numBins, size_t iterations) {
        // Prepare test data
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(-M_PI, M_PI);
        
        std::unique_ptr<double[]> phases(new double[numBins]);
        std::unique_ptr<double[]> workingCopy(new double[numBins]);
        
        for (size_t i = 0; i < numBins; ++i) {
            phases[i] = dist(rng);
        }
        
        // Benchmark scalar implementation
        auto startScalar = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < iterations; ++iter) {
            std::memcpy(workingCopy.get(), phases.get(), numBins * sizeof(double));
            SIMD::Scalar::applyPhaseSmoothing(workingCopy.get(), numBins);
        }
        auto endScalar = std::chrono::high_resolution_clock::now();
        
        // Benchmark SIMD implementation (currently same as scalar)
        auto startSIMD = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < iterations; ++iter) {
            std::memcpy(workingCopy.get(), phases.get(), numBins * sizeof(double));
            SIMD::applyPhaseSmoothingSIMD(workingCopy.get(), numBins);
        }
        auto endSIMD = std::chrono::high_resolution_clock::now();
        
        // Calculate times
        double scalarTime = std::chrono::duration<double, std::milli>(
            endScalar - startScalar).count();
        double simdTime = std::chrono::duration<double, std::milli>(
            endSIMD - startSIMD).count();
        
        return {
            "Phase Smoothing (" + std::to_string(numBins) + " bins)",
            scalarTime,
            simdTime,
            scalarTime / simdTime,
            iterations
        };
    }
    
    static void printResults(const std::vector<BenchmarkResult>& results) {
        const int nameWidth = 60;
        const int numWidth = 12;
        
        std::cout << std::left << std::setw(nameWidth) << "Benchmark"
                  << std::right << std::setw(numWidth) << "Scalar (ms)"
                  << std::setw(numWidth) << "SIMD (ms)"
                  << std::setw(numWidth) << "Speedup"
                  << std::setw(numWidth) << "Iterations"
                  << std::endl;
        
        std::cout << std::string(nameWidth + numWidth * 4, '-') << std::endl;
        
        for (const auto& result : results) {
            std::cout << std::left << std::setw(nameWidth) << result.name
                      << std::right << std::setw(numWidth) << std::fixed 
                      << std::setprecision(2) << result.scalarTime
                      << std::setw(numWidth) << result.simdTime
                      << std::setw(numWidth) << std::setprecision(2) 
                      << result.speedup << "x"
                      << std::setw(numWidth) << result.iterations
                      << std::endl;
        }
        
        std::cout << std::endl;
        
        // Calculate average speedup
        double totalSpeedup = 0.0;
        for (const auto& result : results) {
            totalSpeedup += result.speedup;
        }
        double avgSpeedup = totalSpeedup / results.size();
        
        std::cout << "Average speedup: " << std::fixed << std::setprecision(2) 
                  << avgSpeedup << "x" << std::endl;
    }
};

} // namespace Explorer
} // namespace Acorex

// Standalone function to run benchmarks
void runSIMDBenchmarks() {
    Acorex::Explorer::SIMDBenchmark::runAllBenchmarks();
}