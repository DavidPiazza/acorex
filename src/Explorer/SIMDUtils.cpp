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
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
    // ARM processors support NEON SIMD
    features = features | CpuFeatures::NEON;
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

// ARM NEON implementation of weighted geometric mean
#if defined(__ARM_NEON) || defined(__aarch64__)
void computeWeightedGeometricMeanNEON(
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result) {
    
    const double epsilon = 1e-10;
    const float epsilon_f = 1e-10f;
    
    // Check for single source dominance
    for (size_t i = 0; i < numSources; ++i) {
        if (std::abs(weights[i] - 1.0) < 1e-6) {
            std::memcpy(result, magnitudes[i], numBins * sizeof(double));
            return;
        }
    }
    
    // Process 4 bins at a time using float32x4_t
    // Note: We convert double to float for NEON processing, then back to double
    size_t simdBins = (numBins / 4) * 4;
    
    for (size_t bin = 0; bin < simdBins; bin += 4) {
        float32x4_t logSum = vdupq_n_f32(0.0f);
        
        for (size_t src = 0; src < numSources; ++src) {
            if (weights[src] > epsilon) {
                float weight_f = static_cast<float>(weights[src]);
                float32x4_t weight = vdupq_n_f32(weight_f);
                
                // Convert magnitudes to float
                float mags_f[4];
                for (int i = 0; i < 4; ++i) {
                    mags_f[i] = static_cast<float>(magnitudes[src][bin + i]);
                }
                float32x4_t mag = vld1q_f32(mags_f);
                
                // Add epsilon for numerical stability
                float32x4_t eps = vdupq_n_f32(epsilon_f * weight_f);
                mag = vaddq_f32(mag, eps);
                
                // Compute weighted log using NEON approximation
                float32x4_t logMag = neon_log_ps(mag);
                logSum = vfmaq_f32(logSum, weight, logMag);
            }
        }
        
        // Convert back from log domain
        float32x4_t res = neon_exp_ps(logSum);
        
        // Convert back to double and store
        float res_f[4];
        vst1q_f32(res_f, res);
        for (int i = 0; i < 4; ++i) {
            result[bin + i] = static_cast<double>(res_f[i]);
            if (!std::isfinite(result[bin + i])) {
                result[bin + i] = 0.0;
            }
        }
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

// ARM NEON implementation of weighted circular mean
void computeWeightedCircularMeanNEON(
    const double* const* phases,
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result) {
    
    const double epsilon = 1e-10;
    const float epsilon_f = 1e-10f;
    
    // Process 4 bins at a time
    size_t simdBins = (numBins / 4) * 4;
    
    for (size_t bin = 0; bin < simdBins; bin += 4) {
        float32x4_t sumX = vdupq_n_f32(0.0f);
        float32x4_t sumY = vdupq_n_f32(0.0f);
        float32x4_t totalWeight = vdupq_n_f32(0.0f);
        
        for (size_t src = 0; src < numSources; ++src) {
            if (weights[src] > epsilon) {
                float weight_f = static_cast<float>(weights[src]);
                float32x4_t weight = vdupq_n_f32(weight_f);
                
                // Convert magnitudes and phases to float
                float mags_f[4], phases_f[4];
                for (int i = 0; i < 4; ++i) {
                    mags_f[i] = static_cast<float>(magnitudes[src][bin + i]);
                    phases_f[i] = static_cast<float>(phases[src][bin + i]);
                }
                float32x4_t mag = vld1q_f32(mags_f);
                float32x4_t phase = vld1q_f32(phases_f);
                
                // Combined weight: barycentric weight * magnitude
                float32x4_t w = vmulq_f32(weight, mag);
                
                // Only process if weight is significant
                uint32x4_t mask = vcgtq_f32(w, vdupq_n_f32(epsilon_f));
                w = vbslq_f32(mask, w, vdupq_n_f32(0.0f));
                
                // Compute weighted unit vectors using NEON trig approximations
                float32x4_t cosPhase = neon_cos_ps(phase);
                float32x4_t sinPhase = neon_sin_ps(phase);
                
                sumX = vfmaq_f32(sumX, w, cosPhase);
                sumY = vfmaq_f32(sumY, w, sinPhase);
                totalWeight = vaddq_f32(totalWeight, w);
            }
        }
        
        // Compute mean resultant vector
        float32x4_t invWeight = vdivq_f32(vdupq_n_f32(1.0f), 
                                         vmaxq_f32(totalWeight, vdupq_n_f32(epsilon_f)));
        float32x4_t meanX = vmulq_f32(sumX, invWeight);
        float32x4_t meanY = vmulq_f32(sumY, invWeight);
        
        // Compute circular mean phase using NEON atan2
        float32x4_t meanPhase = neon_atan2_ps(meanY, meanX);
        
        // Convert back to double and store
        float phase_f[4];
        vst1q_f32(phase_f, meanPhase);
        for (int i = 0; i < 4; ++i) {
            result[bin + i] = static_cast<double>(phase_f[i]);
        }
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

// ARM NEON implementation of phase smoothing
void applyPhaseSmoothingNEON(double* phases, size_t numBins) {
    if (numBins < 3) return;
    
    const float32x4_t pi = vdupq_n_f32(M_PI);
    const float32x4_t two_pi = vdupq_n_f32(2.0f * M_PI);
    const float32x4_t neg_two_pi = vdupq_n_f32(-2.0f * M_PI);
    
    // Process 4 bins at a time
    size_t simdBins = ((numBins - 2) / 4) * 4 + 1;
    
    for (size_t i = 1; i < simdBins && i < numBins - 1; i += 4) {
        size_t count = std::min(size_t(4), numBins - 1 - i);
        
        // Load current phases and convert to float
        float curr_f[4], prev_f[4], next_f[4];
        for (size_t j = 0; j < count; ++j) {
            curr_f[j] = static_cast<float>(phases[i + j]);
            prev_f[j] = static_cast<float>(phases[i + j - 1]);
            next_f[j] = static_cast<float>(phases[i + j + 1]);
        }
        
        float32x4_t curr = vld1q_f32(curr_f);
        float32x4_t prev = vld1q_f32(prev_f);
        float32x4_t next = vld1q_f32(next_f);
        
        // Compute differences
        float32x4_t diff_prev = vsubq_f32(curr, prev);
        float32x4_t diff_next = vsubq_f32(next, curr);
        
        // Wrap differences to [-pi, pi]
        uint32x4_t mask_gt_pi = vcgtq_f32(diff_prev, pi);
        uint32x4_t mask_lt_neg_pi = vcltq_f32(diff_prev, vnegq_f32(pi));
        diff_prev = vbslq_f32(mask_gt_pi, vsubq_f32(diff_prev, two_pi), diff_prev);
        diff_prev = vbslq_f32(mask_lt_neg_pi, vaddq_f32(diff_prev, two_pi), diff_prev);
        
        mask_gt_pi = vcgtq_f32(diff_next, pi);
        mask_lt_neg_pi = vcltq_f32(diff_next, vnegq_f32(pi));
        diff_next = vbslq_f32(mask_gt_pi, vsubq_f32(diff_next, two_pi), diff_next);
        diff_next = vbslq_f32(mask_lt_neg_pi, vaddq_f32(diff_next, two_pi), diff_next);
        
        // Apply smoothing if both differences have the same sign
        float32x4_t sign_prev = vbslq_f32(vcltq_f32(diff_prev, vdupq_n_f32(0.0f)), 
                                          vdupq_n_f32(-1.0f), vdupq_n_f32(1.0f));
        float32x4_t sign_next = vbslq_f32(vcltq_f32(diff_next, vdupq_n_f32(0.0f)), 
                                          vdupq_n_f32(-1.0f), vdupq_n_f32(1.0f));
        
        uint32x4_t same_sign = vceqq_f32(sign_prev, sign_next);
        
        // Smooth phase: (prev + curr + next) / 3
        float32x4_t sum = vaddq_f32(vaddq_f32(prev, curr), next);
        float32x4_t smoothed = vmulq_f32(sum, vdupq_n_f32(1.0f / 3.0f));
        
        // Apply smoothing only where signs match
        curr = vbslq_f32(same_sign, smoothed, curr);
        
        // Convert back to double and store
        float result_f[4];
        vst1q_f32(result_f, curr);
        for (size_t j = 0; j < count; ++j) {
            phases[i + j] = static_cast<double>(result_f[j]);
        }
    }
    
    // Handle remaining bins with scalar code
    for (size_t i = simdBins; i < numBins - 1; ++i) {
        double prev = phases[i - 1];
        double curr = phases[i];
        double next = phases[i + 1];
        
        // Compute phase differences
        double diff_prev = curr - prev;
        double diff_next = next - curr;
        
        // Wrap differences to [-pi, pi]
        while (diff_prev > M_PI) diff_prev -= 2.0 * M_PI;
        while (diff_prev < -M_PI) diff_prev += 2.0 * M_PI;
        while (diff_next > M_PI) diff_next -= 2.0 * M_PI;
        while (diff_next < -M_PI) diff_next += 2.0 * M_PI;
        
        // Apply smoothing if both differences have the same sign
        if (diff_prev * diff_next > 0) {
            phases[i] = (prev + curr + next) / 3.0;
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

// Forward declarations for NEON implementations
#if defined(__ARM_NEON) || defined(__aarch64__)
void computeWeightedGeometricMeanNEON(
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result);

void computeWeightedCircularMeanNEON(
    const double* const* phases,
    const double* const* magnitudes,
    const double* weights,
    size_t numSources,
    size_t numBins,
    double* result);

void applyPhaseSmoothingNEON(double* phases, size_t numBins);
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
#if defined(__ARM_NEON) || defined(__aarch64__)
    if (hasFeature(features, CpuFeatures::NEON)) {
        weightedGeometricMean = computeWeightedGeometricMeanNEON;
        weightedCircularMean = computeWeightedCircularMeanNEON;
        phaseSmoothing = applyPhaseSmoothingNEON;
    } else
#endif
#if defined(__AVX__)
    if (hasFeature(features, CpuFeatures::AVX)) {
        weightedGeometricMean = computeWeightedGeometricMeanAVX;
        weightedCircularMean = computeWeightedCircularMeanAVX;
        phaseSmoothing = Scalar::applyPhaseSmoothing;
    } else
#endif
#if defined(__SSE2__)
    if (hasFeature(features, CpuFeatures::SSE2)) {
        weightedGeometricMean = computeWeightedGeometricMeanSSE;
        // No SSE2 circular mean yet, use scalar
        weightedCircularMean = Scalar::computeWeightedCircularMean;
        phaseSmoothing = Scalar::applyPhaseSmoothing;
    } else
#endif
    {
        // Fallback to scalar implementation
        weightedGeometricMean = Scalar::computeWeightedGeometricMean;
        weightedCircularMean = Scalar::computeWeightedCircularMean;
        phaseSmoothing = Scalar::applyPhaseSmoothing;
    }
    
    // Initialize new function pointers
#if defined(__ARM_NEON) || defined(__aarch64__)
    if (hasFeature(features, CpuFeatures::NEON)) {
        // TODO: Add NEON implementations for window/overlap/complex functions
        applyWindow = Scalar::applyWindow;
        overlapAdd = Scalar::overlapAdd;
        complexToMagPhase = Scalar::complexToMagPhase;
        magPhaseToComplex = Scalar::magPhaseToComplex;
        phaseUnwrap = Scalar::phaseUnwrap;
        findSpectralPeaks = Scalar::findSpectralPeaks;
    } else
#endif
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