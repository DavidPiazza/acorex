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
#include <cstring>
#include <algorithm>
#include <vector>

#if defined(_MSC_VER)
#include <intrin.h>
#elif (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#include <cpuid.h>
#endif

#if defined(__ARM_NEON) || defined(__aarch64__)
#include <arm_neon.h>
#endif

namespace Acorex {
namespace Explorer {
namespace SIMD {

CpuFeatures detectCpuFeatures() {
    CpuFeatures features = CpuFeatures::None;
    
#if defined(__aarch64__) || defined(__ARM_NEON)
    // ARM processors don't have SSE/AVX, but we can still use scalar optimizations
    // NEON is the ARM SIMD extension, but we're not implementing it here
    return features;
    
#elif defined(_MSC_VER) && (defined(_M_AMD64) || defined(_M_IX86))
    int cpuInfo[4];
    __cpuid(cpuInfo, 1);
    
    if (cpuInfo[3] & (1 << 26)) features = features | CpuFeatures::SSE2;
    if (cpuInfo[2] & (1 << 19)) features = features | CpuFeatures::SSE41;
    if (cpuInfo[2] & (1 << 28)) features = features | CpuFeatures::AVX;
    if (cpuInfo[2] & (1 << 12)) features = features | CpuFeatures::FMA;
    
    __cpuidex(cpuInfo, 7, 0);
    if (cpuInfo[1] & (1 << 5)) features = features | CpuFeatures::AVX2;
    
#elif (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
    unsigned int eax, ebx, ecx, edx;
    
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        if (edx & (1 << 26)) features = features | CpuFeatures::SSE2;
        if (ecx & (1 << 19)) features = features | CpuFeatures::SSE41;
        if (ecx & (1 << 28)) features = features | CpuFeatures::AVX;
        if (ecx & (1 << 12)) features = features | CpuFeatures::FMA;
    }
    
    if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
        if (ebx & (1 << 5)) features = features | CpuFeatures::AVX2;
    }
#endif
    
    return features;
}

// AVX implementation of weighted geometric mean
#if defined(__AVX__)
void computeWeightedGeometricMeanAVX(
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result) {
    
    const double epsilon = 1e-10;
    
    // Check for single source dominance
    for (size_t i = 0; i < numSources; ++i) {
        if (std::abs(weights[i] - 1.0) < 1e-6) {
            std::memcpy(result, magnitudes[i], numBins * sizeof(double));
            return;
        }
    }
    
    // Process 4 bins at a time with AVX
    size_t simdBins = (numBins / 4) * 4;
    
    for (size_t bin = 0; bin < simdBins; bin += 4) {
        __m256d logSum = _mm256_setzero_pd();
        
        for (size_t src = 0; src < numSources; ++src) {
            if (weights[src] > epsilon) {
                __m256d weight = _mm256_set1_pd(weights[src]);
                __m256d mag = _mm256_loadu_pd(&magnitudes[src][bin]);
                __m256d eps = _mm256_set1_pd(epsilon * weights[src]);
                mag = _mm256_add_pd(mag, eps);
                
                // Compute weighted log (using fast approximation)
                __m256d logMag = log_pd(mag);
                logSum = _mm256_fmadd_pd(weight, logMag, logSum);
            }
        }
        
        // Convert back from log domain
        __m256d res = exp_pd(logSum);
        _mm256_storeu_pd(&result[bin], res);
    }
    
    // Handle remaining bins
    for (size_t bin = simdBins; bin < numBins; ++bin) {
        double logSum = 0.0;
        
        for (size_t src = 0; src < numSources; ++src) {
            if (weights[src] > epsilon) {
                logSum += weights[src] * std::log(magnitudes[src][bin] + epsilon * weights[src]);
            }
        }
        
        result[bin] = std::exp(logSum);
        if (!std::isfinite(result[bin])) {
            result[bin] = 0.0;
        }
    }
}
#endif

// SSE2 implementation of weighted geometric mean
#if defined(__SSE2__)
void computeWeightedGeometricMeanSSE(
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result) {
    
    const double epsilon = 1e-10;
    
    // Check for single source dominance
    for (size_t i = 0; i < numSources; ++i) {
        if (std::abs(weights[i] - 1.0) < 1e-6) {
            std::memcpy(result, magnitudes[i], numBins * sizeof(double));
            return;
        }
    }
    
    // Process 2 bins at a time with SSE2
    size_t simdBins = (numBins / 2) * 2;
    
    for (size_t bin = 0; bin < simdBins; bin += 2) {
        __m128d logSum = _mm_setzero_pd();
        
        for (size_t src = 0; src < numSources; ++src) {
            if (weights[src] > epsilon) {
                __m128d weight = _mm_set1_pd(weights[src]);
                __m128d mag = _mm_loadu_pd(&magnitudes[src][bin]);
                __m128d eps = _mm_set1_pd(epsilon * weights[src]);
                mag = _mm_add_pd(mag, eps);
                
                // Compute weighted log (using scalar fallback for SSE2)
                double m0 = _mm_cvtsd_f64(mag);
                double m1 = _mm_cvtsd_f64(_mm_unpackhi_pd(mag, mag));
                __m128d logMag = _mm_set_pd(std::log(m1), std::log(m0));
                
                logSum = _mm_add_pd(logSum, _mm_mul_pd(weight, logMag));
            }
        }
        
        // Convert back from log domain
        double l0 = _mm_cvtsd_f64(logSum);
        double l1 = _mm_cvtsd_f64(_mm_unpackhi_pd(logSum, logSum));
        __m128d res = _mm_set_pd(std::exp(l1), std::exp(l0));
        _mm_storeu_pd(&result[bin], res);
    }
    
    // Handle remaining bins
    for (size_t bin = simdBins; bin < numBins; ++bin) {
        double logSum = 0.0;
        
        for (size_t src = 0; src < numSources; ++src) {
            if (weights[src] > epsilon) {
                logSum += weights[src] * std::log(magnitudes[src][bin] + epsilon * weights[src]);
            }
        }
        
        result[bin] = std::exp(logSum);
        if (!std::isfinite(result[bin])) {
            result[bin] = 0.0;
        }
    }
}
#endif

// AVX implementation of weighted circular mean
#if defined(__AVX__)
void computeWeightedCircularMeanAVX(
    const double* const* phases,
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result) {
    
    const double epsilon = 1e-10;
    const __m256d eps = _mm256_set1_pd(epsilon);
    const __m256d pi = _mm256_set1_pd(M_PI);
    const __m256d two_pi = _mm256_set1_pd(2.0 * M_PI);
    
    // Process 4 bins at a time
    size_t simdBins = (numBins / 4) * 4;
    
    for (size_t bin = 0; bin < simdBins; bin += 4) {
        __m256d sumX = _mm256_setzero_pd();
        __m256d sumY = _mm256_setzero_pd();
        __m256d totalWeight = _mm256_setzero_pd();
        
        for (size_t src = 0; src < numSources; ++src) {
            if (weights[src] > epsilon) {
                __m256d weight = _mm256_set1_pd(weights[src]);
                __m256d mag = _mm256_loadu_pd(&magnitudes[src][bin]);
                __m256d phase = _mm256_loadu_pd(&phases[src][bin]);
                
                // Combined weight: barycentric weight * magnitude
                __m256d w = _mm256_mul_pd(weight, mag);
                
                // Only process if weight is significant
                __m256d mask = _mm256_cmp_pd(w, eps, _CMP_GT_OQ);
                w = _mm256_and_pd(w, mask);
                
                // Compute weighted unit vectors
                // Note: Using scalar sin/cos for now - could be optimized with Taylor series
                alignas(32) double w_array[4];
                alignas(32) double phase_array[4];
                _mm256_store_pd(w_array, w);
                _mm256_store_pd(phase_array, phase);
                
                __m256d cosPhase = _mm256_set_pd(
                    std::cos(phase_array[3]), std::cos(phase_array[2]),
                    std::cos(phase_array[1]), std::cos(phase_array[0])
                );
                __m256d sinPhase = _mm256_set_pd(
                    std::sin(phase_array[3]), std::sin(phase_array[2]),
                    std::sin(phase_array[1]), std::sin(phase_array[0])
                );
                
                sumX = _mm256_fmadd_pd(w, cosPhase, sumX);
                sumY = _mm256_fmadd_pd(w, sinPhase, sumY);
                totalWeight = _mm256_add_pd(totalWeight, w);
            }
        }
        
        // Compute mean resultant vector
        __m256d mask = _mm256_cmp_pd(totalWeight, eps, _CMP_GT_OQ);
        __m256d invWeight = _mm256_div_pd(_mm256_set1_pd(1.0), 
                                         _mm256_max_pd(totalWeight, eps));
        __m256d meanX = _mm256_mul_pd(sumX, invWeight);
        __m256d meanY = _mm256_mul_pd(sumY, invWeight);
        
        // Compute circular mean phase
        alignas(32) double meanX_array[4];
        alignas(32) double meanY_array[4];
        _mm256_store_pd(meanX_array, meanX);
        _mm256_store_pd(meanY_array, meanY);
        
        __m256d meanPhase = _mm256_set_pd(
            std::atan2(meanY_array[3], meanX_array[3]),
            std::atan2(meanY_array[2], meanX_array[2]),
            std::atan2(meanY_array[1], meanX_array[1]),
            std::atan2(meanY_array[0], meanX_array[0])
        );
        
        _mm256_storeu_pd(&result[bin], meanPhase);
    }
    
    // Handle remaining bins with scalar code
    for (size_t bin = simdBins; bin < numBins; ++bin) {
        double sumX = 0.0;
        double sumY = 0.0;
        double totalWeight = 0.0;
        
        for (size_t src = 0; src < numSources; ++src) {
            if (weights[src] > epsilon && magnitudes[src][bin] > epsilon) {
                double w = weights[src] * magnitudes[src][bin];
                double phase = phases[src][bin];
                
                sumX += w * std::cos(phase);
                sumY += w * std::sin(phase);
                totalWeight += w;
            }
        }
        
        if (totalWeight > epsilon) {
            double meanX = sumX / totalWeight;
            double meanY = sumY / totalWeight;
            result[bin] = std::atan2(meanY, meanX);
        } else {
            result[bin] = 0.0;
        }
    }
}

// AVX implementation of window application
void applyWindowAVX(const double* input, const double* window, size_t windowSize, double* output) {
    size_t simdSize = (windowSize / 4) * 4;
    
    // Process 4 samples at a time
    for (size_t i = 0; i < simdSize; i += 4) {
        __m256d in = _mm256_loadu_pd(&input[i]);
        __m256d win = _mm256_loadu_pd(&window[i]);
        __m256d result = _mm256_mul_pd(in, win);
        _mm256_storeu_pd(&output[i], result);
    }
    
    // Handle remaining samples
    for (size_t i = simdSize; i < windowSize; ++i) {
        output[i] = input[i] * window[i];
    }
}

// AVX implementation of overlap-add
void overlapAddAVX(const double* input, double* accumulator, size_t frameSize, size_t hopSize) {
    size_t simdSize = (frameSize / 4) * 4;
    
    // Process 4 samples at a time
    for (size_t i = 0; i < simdSize; i += 4) {
        __m256d in = _mm256_loadu_pd(&input[i]);
        __m256d acc = _mm256_loadu_pd(&accumulator[i]);
        __m256d result = _mm256_add_pd(acc, in);
        _mm256_storeu_pd(&accumulator[i], result);
    }
    
    // Handle remaining samples
    for (size_t i = simdSize; i < frameSize; ++i) {
        accumulator[i] += input[i];
    }
}

// AVX implementation of complex to magnitude/phase conversion
void complexToMagPhaseAVX(const std::complex<double>* spectrum, size_t numBins,
                         double* magnitudes, double* phases) {
    // Process 2 complex numbers at a time (4 doubles)
    size_t simdBins = (numBins / 2) * 2;
    
    for (size_t i = 0; i < simdBins; i += 2) {
        // Load 2 complex numbers (real0, imag0, real1, imag1)
        __m256d complexVec = _mm256_loadu_pd(reinterpret_cast<const double*>(&spectrum[i]));
        
        // Extract real and imaginary parts
        __m256d realParts = _mm256_shuffle_pd(complexVec, complexVec, 0x0); // real0, real0, real1, real1
        __m256d imagParts = _mm256_shuffle_pd(complexVec, complexVec, 0xF); // imag0, imag0, imag1, imag1
        
        // Compute magnitudes: sqrt(real^2 + imag^2)
        __m256d real2 = _mm256_mul_pd(realParts, realParts);
        __m256d imag2 = _mm256_mul_pd(imagParts, imagParts);
        __m256d sum = _mm256_add_pd(real2, imag2);
        __m256d mags = _mm256_sqrt_pd(sum);
        
        // Store magnitudes (we only need elements 0 and 2)
        magnitudes[i] = mags[0];
        magnitudes[i+1] = mags[2];
        
        // Compute phases using atan2 (scalar fallback for now)
        phases[i] = std::atan2(spectrum[i].imag(), spectrum[i].real());
        phases[i+1] = std::atan2(spectrum[i+1].imag(), spectrum[i+1].real());
    }
    
    // Handle remaining bins
    for (size_t i = simdBins; i < numBins; ++i) {
        magnitudes[i] = std::abs(spectrum[i]);
        phases[i] = std::arg(spectrum[i]);
    }
}
#endif

// Scalar implementations
namespace Scalar {

void computeWeightedGeometricMean(
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result) {
    
    const double epsilon = 1e-10;
    
    // Check for single source dominance
    for (size_t i = 0; i < numSources; ++i) {
        if (std::abs(weights[i] - 1.0) < 1e-6) {
            std::memcpy(result, magnitudes[i], numBins * sizeof(double));
            return;
        }
    }
    
    // Compute weighted sum in log domain
    for (size_t bin = 0; bin < numBins; ++bin) {
        double logSum = 0.0;
        
        for (size_t src = 0; src < numSources; ++src) {
            if (weights[src] > epsilon) {
                logSum += weights[src] * std::log(magnitudes[src][bin] + epsilon * weights[src]);
            }
        }
        
        result[bin] = std::exp(logSum);
        if (!std::isfinite(result[bin])) {
            result[bin] = 0.0;
        }
    }
}

void computeWeightedCircularMean(
    const double* const* phases,
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result) {
    
    const double epsilon = 1e-10;
    
    for (size_t bin = 0; bin < numBins; ++bin) {
        double sumX = 0.0;
        double sumY = 0.0;
        double totalWeight = 0.0;
        
        for (size_t src = 0; src < numSources; ++src) {
            if (weights[src] > epsilon && magnitudes[src][bin] > epsilon) {
                double w = weights[src] * magnitudes[src][bin];
                double phase = phases[src][bin];
                
                sumX += w * std::cos(phase);
                sumY += w * std::sin(phase);
                totalWeight += w;
            }
        }
        
        if (totalWeight > epsilon) {
            double meanX = sumX / totalWeight;
            double meanY = sumY / totalWeight;
            result[bin] = std::atan2(meanY, meanX);
            
            // Handle phase unwrapping
            if (bin > 0) {
                double prevPhase = result[bin - 1];
                double phaseDiff = result[bin] - prevPhase;
                
                while (phaseDiff > M_PI) {
                    result[bin] -= 2 * M_PI;
                    phaseDiff = result[bin] - prevPhase;
                }
                while (phaseDiff < -M_PI) {
                    result[bin] += 2 * M_PI;
                    phaseDiff = result[bin] - prevPhase;
                }
            }
        } else {
            result[bin] = 0.0;
        }
    }
}

void applyPhaseSmoothing(double* phases, size_t numBins) {
    if (numBins < 3) return;
    
    const double maxPhaseJump = M_PI / 2;
    const int smoothingWindow = 5;
    
    // Detect discontinuities
    std::vector<bool> discontinuities(numBins, false);
    
    for (size_t i = 0; i < numBins - 1; ++i) {
        double diff = phases[i + 1] - phases[i];
        
        while (diff > M_PI) diff -= 2 * M_PI;
        while (diff < -M_PI) diff += 2 * M_PI;
        
        if (std::abs(diff) > maxPhaseJump) {
            discontinuities[i] = true;
            discontinuities[i + 1] = true;
        }
    }
    
    // Apply median filtering near discontinuities
    std::vector<double> smoothedPhases(phases, phases + numBins);
    
    for (size_t i = 0; i < numBins; ++i) {
        if (discontinuities[i]) {
            std::vector<double> window;
            
            for (int w = -smoothingWindow/2; w <= smoothingWindow/2; ++w) {
                int idx = static_cast<int>(i) + w;
                
                if (idx < 0) idx = -idx;
                if (idx >= static_cast<int>(numBins)) {
                    idx = 2 * static_cast<int>(numBins) - idx - 2;
                }
                
                if (idx >= 0 && idx < static_cast<int>(numBins)) {
                    double phase = phases[idx];
                    double centerPhase = phases[i];
                    double diff = phase - centerPhase;
                    
                    while (diff > M_PI) diff -= 2 * M_PI;
                    while (diff < -M_PI) diff += 2 * M_PI;
                    
                    window.push_back(centerPhase + diff);
                }
            }
            
            if (!window.empty()) {
                std::sort(window.begin(), window.end());
                smoothedPhases[i] = window[window.size() / 2];
            }
        }
    }
    
    // Ensure phase continuity
    for (size_t i = 1; i < numBins; ++i) {
        double diff = smoothedPhases[i] - smoothedPhases[i - 1];
        
        while (diff > M_PI) {
            smoothedPhases[i] -= 2 * M_PI;
            diff = smoothedPhases[i] - smoothedPhases[i - 1];
        }
        while (diff < -M_PI) {
            smoothedPhases[i] += 2 * M_PI;
            diff = smoothedPhases[i] - smoothedPhases[i - 1];
        }
    }
    
    std::memcpy(phases, smoothedPhases.data(), numBins * sizeof(double));
}

void applyWindow(const double* input, const double* window, size_t windowSize, double* output) {
    for (size_t i = 0; i < windowSize; ++i) {
        output[i] = input[i] * window[i];
    }
}

void overlapAdd(const double* input, double* accumulator, size_t frameSize, size_t hopSize) {
    for (size_t i = 0; i < frameSize; ++i) {
        accumulator[i] += input[i];
    }
}

void complexToMagPhase(const std::complex<double>* spectrum, size_t numBins,
                      double* magnitudes, double* phases) {
    for (size_t i = 0; i < numBins; ++i) {
        magnitudes[i] = std::abs(spectrum[i]);
        phases[i] = std::arg(spectrum[i]);
    }
}

void magPhaseToComplex(const double* magnitudes, const double* phases,
                      size_t numBins, std::complex<double>* spectrum) {
    for (size_t i = 0; i < numBins; ++i) {
        spectrum[i] = std::polar(magnitudes[i], phases[i]);
    }
}

void phaseUnwrap(double* phases, size_t numBins) {
    if (numBins < 2) return;
    
    for (size_t i = 1; i < numBins; ++i) {
        double diff = phases[i] - phases[i-1];
        
        // Unwrap phase jumps greater than π
        while (diff > M_PI) {
            phases[i] -= 2.0 * M_PI;
            diff = phases[i] - phases[i-1];
        }
        while (diff < -M_PI) {
            phases[i] += 2.0 * M_PI;
            diff = phases[i] - phases[i-1];
        }
    }
}

void findSpectralPeaks(const double* magnitudes, size_t numBins, double threshold,
                      size_t* peaks, size_t& numPeaks) {
    numPeaks = 0;
    const size_t maxPeaks = 1024; // Reasonable limit
    
    // Skip DC and Nyquist bins
    for (size_t i = 1; i < numBins - 1; ++i) {
        // Check if this bin is a local maximum above threshold
        if (magnitudes[i] > threshold &&
            magnitudes[i] > magnitudes[i-1] &&
            magnitudes[i] > magnitudes[i+1]) {
            
            peaks[numPeaks++] = i;
            
            if (numPeaks >= maxPeaks) break;
        }
    }
}

} // namespace Scalar

// Runtime dispatcher implementation
SIMDDispatcher::SIMDDispatcher() {
    features = detectCpuFeatures();
    
    // Select best available implementation
#if defined(__AVX__)
    if (hasFeature(features, CpuFeatures::AVX)) {
        weightedGeometricMean = computeWeightedGeometricMeanAVX;
        weightedCircularMean = computeWeightedCircularMeanAVX;
    } else
#endif
#if defined(__SSE2__)
    if (hasFeature(features, CpuFeatures::SSE2)) {
        weightedGeometricMean = computeWeightedGeometricMeanSSE;
        // No SSE2 circular mean yet, use scalar
        weightedCircularMean = Scalar::computeWeightedCircularMean;
    } else
#endif
    {
        // Fallback to scalar implementation
        weightedGeometricMean = Scalar::computeWeightedGeometricMean;
        weightedCircularMean = Scalar::computeWeightedCircularMean;
    }
    
    // Phase smoothing is currently scalar only
    phaseSmoothing = Scalar::applyPhaseSmoothing;
    
    // Initialize new function pointers
#if defined(__AVX__)
    if (hasFeature(features, CpuFeatures::AVX)) {
        applyWindow = applyWindowAVX;
        overlapAdd = overlapAddAVX;
        complexToMagPhase = complexToMagPhaseAVX;
        // TODO: Add AVX implementations for remaining functions
        magPhaseToComplex = Scalar::magPhaseToComplex;
        phaseUnwrap = Scalar::phaseUnwrap;
        findSpectralPeaks = Scalar::findSpectralPeaks;
    } else
#endif
    {
        // Fallback to scalar implementations
        applyWindow = Scalar::applyWindow;
        overlapAdd = Scalar::overlapAdd;
        complexToMagPhase = Scalar::complexToMagPhase;
        magPhaseToComplex = Scalar::magPhaseToComplex;
        phaseUnwrap = Scalar::phaseUnwrap;
        findSpectralPeaks = Scalar::findSpectralPeaks;
    }
}

// Public API functions that use runtime dispatch
void computeWeightedGeometricMeanSIMD(
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result) {
    
    const auto& dispatcher = SIMDDispatcher::getInstance();
    dispatcher.weightedGeometricMean(magnitudes, weights, numSources, numBins, result);
}

void computeWeightedCircularMeanSIMD(
    const double* const* phases,
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result) {
    
    const auto& dispatcher = SIMDDispatcher::getInstance();
    dispatcher.weightedCircularMean(phases, magnitudes, weights, numSources, numBins, result);
}

void applyPhaseSmoothingSIMD(double* phases, size_t numBins) {
    const auto& dispatcher = SIMDDispatcher::getInstance();
    dispatcher.phaseSmoothing(phases, numBins);
}

void applyWindowSIMD(const double* input, const double* window, size_t windowSize, double* output) {
    const auto& dispatcher = SIMDDispatcher::getInstance();
    dispatcher.applyWindow(input, window, windowSize, output);
}

void overlapAddSIMD(const double* input, double* accumulator, size_t frameSize, size_t hopSize) {
    const auto& dispatcher = SIMDDispatcher::getInstance();
    dispatcher.overlapAdd(input, accumulator, frameSize, hopSize);
}

void complexToMagPhaseSIMD(const std::complex<double>* spectrum, size_t numBins, 
                          double* magnitudes, double* phases) {
    const auto& dispatcher = SIMDDispatcher::getInstance();
    dispatcher.complexToMagPhase(spectrum, numBins, magnitudes, phases);
}

void magPhaseToComplexSIMD(const double* magnitudes, const double* phases, 
                          size_t numBins, std::complex<double>* spectrum) {
    const auto& dispatcher = SIMDDispatcher::getInstance();
    dispatcher.magPhaseToComplex(magnitudes, phases, numBins, spectrum);
}

void phaseUnwrapSIMD(double* phases, size_t numBins) {
    const auto& dispatcher = SIMDDispatcher::getInstance();
    dispatcher.phaseUnwrap(phases, numBins);
}

void findSpectralPeaksSIMD(const double* magnitudes, size_t numBins, double threshold,
                          size_t* peaks, size_t& numPeaks) {
    const auto& dispatcher = SIMDDispatcher::getInstance();
    dispatcher.findSpectralPeaks(magnitudes, numBins, threshold, peaks, numPeaks);
}

} // namespace SIMD
} // namespace Explorer
} // namespace Acorex