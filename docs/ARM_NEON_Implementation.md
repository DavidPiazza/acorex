# ARM NEON SIMD Implementation

## Overview

This document describes the ARM NEON SIMD optimizations implemented for the ACorEx audio processing system. These optimizations provide significant performance improvements on ARM-based platforms including Apple Silicon (M1/M2/M3) and other ARM64 processors.

## Implementation Details

### 1. CPU Feature Detection
- Added `NEON` to the `CpuFeatures` enum in `SIMDUtils.hpp`
- Updated `detectCpuFeatures()` to detect ARM NEON support on ARM/ARM64 platforms
- NEON is automatically detected on all ARM64 platforms and ARM platforms with NEON support

### 2. Memory Alignment
- Added `ARM_NEON_ALIGN` macro for 16-byte alignment requirements
- NEON requires 16-byte alignment for optimal performance with 128-bit vectors

### 3. Polynomial Approximations for Transcendental Functions

Implemented fast NEON approximations for:
- **`neon_exp_ps()`**: Exponential function using range reduction and polynomial approximation
- **`neon_log_ps()`**: Logarithm using mantissa/exponent decomposition
- **`neon_sin_ps()`**: Sine function with range reduction to [-π, π] and Taylor series
- **`neon_cos_ps()`**: Cosine function with range reduction and Taylor series
- **`neon_atan2_ps()`**: Two-argument arctangent using polynomial approximation

### 4. Optimized Audio Processing Functions

#### Weighted Geometric Mean (`computeWeightedGeometricMeanNEON`)
- Processes 4 frequency bins simultaneously using float32x4_t vectors
- Converts double precision to single precision for NEON processing
- Uses NEON log/exp approximations for performance
- Falls back to scalar code for remaining bins

#### Weighted Circular Mean (`computeWeightedCircularMeanNEON`)
- Computes weighted average of phase values on unit circle
- Uses NEON sin/cos approximations for vector decomposition
- Implements NEON atan2 for phase reconstruction
- Handles phase wrapping and numerical stability

#### Phase Smoothing (`applyPhaseSmoothingNEON`)
- Smooths phase discontinuities in frequency domain
- Processes 4 bins at a time with NEON vectorization
- Detects and corrects phase jumps > π
- Applies smoothing only when adjacent phase differences have same sign

### 5. Runtime Dispatcher Integration
- Updated `SIMDDispatcher` constructor to select NEON implementations when available
- NEON functions have priority over SSE/AVX on ARM platforms
- Automatic fallback to scalar implementations if NEON not available

## Performance Expectations

Based on the implementation, expected performance improvements on ARM platforms:
- **Weighted Geometric Mean**: 2-3x speedup
- **Weighted Circular Mean**: 2-4x speedup  
- **Phase Smoothing**: 2-3x speedup

Actual performance gains depend on:
- ARM processor generation (ARMv7 with NEON vs ARMv8/ARM64)
- Memory bandwidth and cache configuration
- Data size and alignment
- Compiler optimization settings

## Testing

### Unit Test
Run the NEON-specific test:
```bash
make test-neon
```

This test:
- Detects CPU features
- Benchmarks NEON vs scalar implementations
- Verifies numerical accuracy
- Reports speedup factors

### Full Benchmark
Run comprehensive SIMD benchmarks:
```bash
make benchmark-simd
```

## Build Configuration

The NEON optimizations are automatically enabled when building on ARM platforms. No special configuration is required.

### Compiler Flags
- ARM32: Requires `-mfpu=neon` (automatically detected)
- ARM64: NEON is always available
- All platforms: Use `-O3` for best performance

## Limitations and Future Work

### Current Limitations
1. Uses single precision (float32) internally, which may introduce small numerical errors
2. Some functions (window, overlap-add, complex operations) still use scalar fallbacks
3. Transcendental function approximations trade accuracy for speed

### Future Improvements
1. Implement NEON versions of remaining functions
2. Add float64x2_t support for ARM64 (double precision NEON)
3. Optimize memory access patterns for better cache utilization
4. Add NEON-specific prefetch instructions
5. Implement fused operations for common patterns

## Compatibility

The implementation is compatible with:
- ARMv7 with NEON extension
- ARMv8/ARM64 (all versions)
- Apple Silicon (M1, M2, M3)
- Raspberry Pi 3/4
- Android ARM devices

The code automatically falls back to scalar implementations on:
- Older ARM processors without NEON
- x86/x64 platforms (uses SSE/AVX instead)