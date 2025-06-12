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

#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <complex>

#if defined(__SSE2__)
#include <emmintrin.h>
#endif

#if defined(__SSE4_1__)
#include <smmintrin.h>
#endif

#if defined(__AVX__)
#include <immintrin.h>
#endif

#if defined(__AVX2__)
#include <immintrin.h>
#endif

#if defined(__FMA__)
#include <immintrin.h>
#endif

namespace Acorex {
namespace Explorer {
namespace SIMD {

enum class CpuFeatures : uint32_t {
    None = 0,
    SSE2 = 1 << 0,
    SSE41 = 1 << 1,
    AVX = 1 << 2,
    AVX2 = 1 << 3,
    FMA = 1 << 4
};

inline CpuFeatures operator|(CpuFeatures a, CpuFeatures b) {
    return static_cast<CpuFeatures>(static_cast<uint32_t>(a) | static_cast<uint32_t>(b));
}

inline CpuFeatures operator&(CpuFeatures a, CpuFeatures b) {
    return static_cast<CpuFeatures>(static_cast<uint32_t>(a) & static_cast<uint32_t>(b));
}

inline bool hasFeature(CpuFeatures features, CpuFeatures feature) {
    return (features & feature) == feature;
}

CpuFeatures detectCpuFeatures();

// SIMD-optimized weighted geometric mean computation
void computeWeightedGeometricMeanSIMD(
    const double* const* magnitudes,  // Array of pointers to magnitude arrays
    const double* weights,            // Barycentric weights
    size_t numSources,               // Number of sources
    size_t numBins,                  // Number of frequency bins
    double* result                   // Output array
);

// SIMD-optimized weighted circular mean computation
void computeWeightedCircularMeanSIMD(
    const double* const* phases,      // Array of pointers to phase arrays
    const double* const* magnitudes,  // Array of pointers to magnitude arrays
    const double* weights,            // Barycentric weights
    size_t numSources,               // Number of sources
    size_t numBins,                  // Number of frequency bins
    double* result                   // Output array
);

// SIMD-optimized phase smoothing
void applyPhaseSmoothingSIMD(
    double* phases,                  // Phase array to smooth (modified in-place)
    size_t numBins                   // Number of frequency bins
);

// SIMD-optimized window application
void applyWindowSIMD(
    const double* input,             // Input audio buffer
    const double* window,            // Window function
    size_t windowSize,               // Size of window
    double* output                   // Output buffer
);

// SIMD-optimized overlap-add accumulation
void overlapAddSIMD(
    const double* input,             // Input windowed frame
    double* accumulator,             // Accumulator buffer
    size_t frameSize,                // Size of frame
    size_t hopSize                   // Hop size
);

// SIMD-optimized complex spectrum to magnitude/phase conversion
void complexToMagPhaseSIMD(
    const std::complex<double>* spectrum,  // Complex spectrum
    size_t numBins,                       // Number of frequency bins
    double* magnitudes,                   // Output magnitudes
    double* phases                        // Output phases
);

// SIMD-optimized magnitude/phase to complex spectrum conversion
void magPhaseToComplexSIMD(
    const double* magnitudes,             // Input magnitudes
    const double* phases,                 // Input phases
    size_t numBins,                      // Number of frequency bins
    std::complex<double>* spectrum       // Output complex spectrum
);

// SIMD-optimized phase unwrapping
void phaseUnwrapSIMD(
    double* phases,                      // Phase array (modified in-place)
    size_t numBins                       // Number of frequency bins
);

// SIMD-optimized spectral peak detection
void findSpectralPeaksSIMD(
    const double* magnitudes,            // Magnitude spectrum
    size_t numBins,                     // Number of frequency bins
    double threshold,                   // Peak threshold
    size_t* peaks,                      // Output peak indices
    size_t& numPeaks                   // Number of peaks found
);

// Scalar fallback implementations
namespace Scalar {
    void computeWeightedGeometricMean(
        const double* const* magnitudes,
        const double* weights,
        size_t numSources,
        size_t numBins,
        double* result
    );
    
    void computeWeightedCircularMean(
        const double* const* phases,
        const double* const* magnitudes,
        const double* weights,
        size_t numSources,
        size_t numBins,
        double* result
    );
    
    void applyPhaseSmoothing(
        double* phases,
        size_t numBins
    );
    
    void applyWindow(
        const double* input,
        const double* window,
        size_t windowSize,
        double* output
    );
    
    void overlapAdd(
        const double* input,
        double* accumulator,
        size_t frameSize,
        size_t hopSize
    );
    
    void complexToMagPhase(
        const std::complex<double>* spectrum,
        size_t numBins,
        double* magnitudes,
        double* phases
    );
    
    void magPhaseToComplex(
        const double* magnitudes,
        const double* phases,
        size_t numBins,
        std::complex<double>* spectrum
    );
    
    void phaseUnwrap(
        double* phases,
        size_t numBins
    );
    
    void findSpectralPeaks(
        const double* magnitudes,
        size_t numBins,
        double threshold,
        size_t* peaks,
        size_t& numPeaks
    );
}

// Function pointer types for runtime dispatch
using WeightedGeometricMeanFunc = void(*)(
    const double* const*, const double*, size_t, size_t, double*);
using WeightedCircularMeanFunc = void(*)(
    const double* const*, const double* const*, const double*, size_t, size_t, double*);
using PhaseSmoothingFunc = void(*)(double*, size_t);
using WindowFunc = void(*)(const double*, const double*, size_t, double*);
using OverlapAddFunc = void(*)(const double*, double*, size_t, size_t);
using ComplexToMagPhaseFunc = void(*)(const std::complex<double>*, size_t, double*, double*);
using MagPhaseToComplexFunc = void(*)(const double*, const double*, size_t, std::complex<double>*);
using PhaseUnwrapFunc = void(*)(double*, size_t);
using SpectralPeaksFunc = void(*)(const double*, size_t, double, size_t*, size_t&);

// Runtime dispatch structure
struct SIMDDispatcher {
    WeightedGeometricMeanFunc weightedGeometricMean;
    WeightedCircularMeanFunc weightedCircularMean;
    PhaseSmoothingFunc phaseSmoothing;
    WindowFunc applyWindow;
    OverlapAddFunc overlapAdd;
    ComplexToMagPhaseFunc complexToMagPhase;
    MagPhaseToComplexFunc magPhaseToComplex;
    PhaseUnwrapFunc phaseUnwrap;
    SpectralPeaksFunc findSpectralPeaks;
    CpuFeatures features;
    
    SIMDDispatcher();
    
    static const SIMDDispatcher& getInstance() {
        static SIMDDispatcher instance;
        return instance;
    }
};

// Helper functions for SIMD operations
#if defined(__AVX__)
inline __m256d exp_pd(__m256d x) {
    // Fast approximation of exp() for AVX
    // Based on the identity: exp(x) = exp(x/log(2)) * 2^(x/log(2))
    const __m256d log2e = _mm256_set1_pd(1.4426950408889634);
    const __m256d ln2 = _mm256_set1_pd(0.6931471805599453);
    
    // Compute x / log(2)
    __m256d a = _mm256_mul_pd(x, log2e);
    
    // Extract integer and fractional parts
    __m256d k = _mm256_round_pd(a, _MM_FROUND_TO_NEAREST_INT);
    __m256d f = _mm256_sub_pd(a, k);
    
    // Compute 2^k
    __m256i ki = _mm256_cvtpd_epi32(k);
    __m256i exp_k = _mm256_slli_epi32(_mm256_add_epi32(ki, _mm256_set1_epi32(1023)), 20);
    __m256d pow2_k = _mm256_castsi256_pd(_mm256_shuffle_epi32(exp_k, _MM_SHUFFLE(1,3,0,2)));
    
    // Polynomial approximation of exp(f * ln(2))
    __m256d p = _mm256_fmadd_pd(f, _mm256_set1_pd(0.0013530747), _mm256_set1_pd(0.0096824656));
    p = _mm256_fmadd_pd(p, f, _mm256_set1_pd(0.0555044437));
    p = _mm256_fmadd_pd(p, f, _mm256_set1_pd(0.2402264689));
    p = _mm256_fmadd_pd(p, f, _mm256_set1_pd(0.6931471806));
    p = _mm256_fmadd_pd(p, f, _mm256_set1_pd(1.0));
    
    return _mm256_mul_pd(pow2_k, p);
}

inline __m256d log_pd(__m256d x) {
    // Fast approximation of log() for AVX
    const __m256d ln2 = _mm256_set1_pd(0.6931471805599453);
    
    // Extract exponent and mantissa
    __m256i xi = _mm256_castpd_si256(x);
    __m256i exp = _mm256_srli_epi64(xi, 52);
    exp = _mm256_sub_epi64(exp, _mm256_set1_epi64x(1023));
    __m256d e = _mm256_cvtepi64_pd(exp);
    
    // Normalize mantissa to [1, 2)
    xi = _mm256_and_si256(xi, _mm256_set1_epi64x(0xFFFFFFFFFFFFF));
    xi = _mm256_or_si256(xi, _mm256_set1_epi64x(0x3FF0000000000000));
    __m256d m = _mm256_castsi256_pd(xi);
    
    // Polynomial approximation of log(m) for m in [1, 2)
    __m256d p = _mm256_set1_pd(-0.6491106);
    p = _mm256_fmadd_pd(p, m, _mm256_set1_pd(3.0734525));
    p = _mm256_fmadd_pd(p, m, _mm256_set1_pd(-5.4195777));
    p = _mm256_fmadd_pd(p, m, _mm256_set1_pd(3.9754935));
    
    // Combine: log(x) = log(m * 2^e) = log(m) + e * log(2)
    return _mm256_fmadd_pd(e, ln2, p);
}
#endif

#if defined(__SSE2__)
inline __m128d exp_pd_sse(__m128d x) {
    // SSE2 version of exp approximation
    const __m128d log2e = _mm_set1_pd(1.4426950408889634);
    const __m128d ln2 = _mm_set1_pd(0.6931471805599453);
    
    __m128d a = _mm_mul_pd(x, log2e);
    __m128d k = _mm_round_pd(a, _MM_FROUND_TO_NEAREST_INT);
    __m128d f = _mm_sub_pd(a, k);
    
    // Convert to integer for bit manipulation
    __m128i ki = _mm_cvtpd_epi32(k);
    ki = _mm_shuffle_epi32(ki, _MM_SHUFFLE(1,1,0,0));
    __m128i exp_k = _mm_slli_epi32(_mm_add_epi32(ki, _mm_set1_epi32(1023)), 20);
    __m128i exp_k_hi = _mm_shuffle_epi32(exp_k, _MM_SHUFFLE(3,3,2,2));
    __m128d pow2_k = _mm_castsi128_pd(_mm_unpacklo_epi32(_mm_setzero_si128(), exp_k_hi));
    
    // Polynomial approximation
    __m128d p = _mm_set1_pd(0.0013530747);
    p = _mm_add_pd(_mm_mul_pd(p, f), _mm_set1_pd(0.0096824656));
    p = _mm_add_pd(_mm_mul_pd(p, f), _mm_set1_pd(0.0555044437));
    p = _mm_add_pd(_mm_mul_pd(p, f), _mm_set1_pd(0.2402264689));
    p = _mm_add_pd(_mm_mul_pd(p, f), _mm_set1_pd(0.6931471806));
    p = _mm_add_pd(_mm_mul_pd(p, f), _mm_set1_pd(1.0));
    
    return _mm_mul_pd(pow2_k, p);
}
#endif

} // namespace SIMD
} // namespace Explorer
} // namespace Acorex