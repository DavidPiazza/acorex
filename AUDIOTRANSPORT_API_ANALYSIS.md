# AudioTransport API Analysis

## Overview

The AudioTransport implementation in this project is based on FluCoMa's AudioTransport algorithm, which implements the paper "Audio Transport: a generalized portamento via optimal transport" by T. Henderson and J. Solomon (DAFX 2019).

## Key Files Found

### 1. Core FluCoMa AudioTransport Implementation
- **Path**: `/deps/flucoma-core/algorithms/public/AudioTransport.hpp`
- **Purpose**: Original pairwise audio morphing using optimal transport
- **Key Components**:
  - `SpectralMass` struct - represents spectral segments with mass
  - `segmentSpectrum()` - segments spectrum into masses based on reassigned frequencies
  - `computeTransportMatrix()` - computes optimal transport plan between two mass distributions
  - `placeMass()` - places interpolated spectral mass at target location
  - `interpolate()` - main morphing function using STFT analysis/synthesis

### 2. N-Way AudioTransport Extension
- **Path**: `/src/Explorer/AudioTransportN.hpp` and `AudioTransportN.cpp`
- **Purpose**: Extends FluCoMa's AudioTransport to support N-way (2-8 sources) morphing
- **Key Features**:
  - `BarycentricWeights` struct for N-way weight management
  - `processFrameN()` - N-way morphing with barycentric coordinates
  - Magnitude interpolation using weighted geometric mean
  - Phase interpolation using circular statistics
  - Fallback to pairwise when N=2 for efficiency

### 3. AudioPlayback Integration
- **Path**: `/src/Explorer/AudioPlayback.cpp` and `AudioPlayback.h`
- **Integration Points**:
  - `mAudioTransport` - standard FluCoMa transport for pairwise morphing
  - `mAudioTransportN` - N-way transport for multi-source morphing
  - `ProcessMorphFrame()` - handles pairwise morphing during crossfades
  - `ProcessNWayMorphFrame()` - handles N-way morphing (partially implemented)
  - Thread-safe morph target and weight management

### 4. FluCoMa Client Wrapper
- **Path**: `/deps/flucoma-core/clients/rt/AudioTransportClient.hpp`
- **Purpose**: Real-time client wrapper for AudioTransport algorithm
- **Features**:
  - Parameter management (interpolation weight, FFT settings)
  - Buffered processing for real-time operation
  - Window normalization

## Algorithm Details

### Optimal Transport Algorithm (from AudioTransport.hpp)

1. **Spectral Segmentation** (`segmentSpectrum`):
   - Uses reassigned frequencies to find spectral peaks/valleys
   - Segments spectrum into "masses" at local minima of reassigned frequency
   - Each mass has: startBin, centerBin, endBin, and total mass (normalized magnitude sum)

2. **Transport Matrix Computation** (`computeTransportMatrix`):
   - Implements 1D Wasserstein distance/Earth Mover's Distance
   - Creates mapping between masses in source A and source B
   - Uses greedy algorithm to match masses in order
   - Output: list of (indexA, indexB, mass) tuples

3. **Interpolation Process**:
   - For each transport tuple:
     - Interpolate bin position between source masses
     - Interpolate frequency using reassigned frequencies
     - Place both source masses at interpolated position with weighted amplitudes
   - Phase tracking for coherent synthesis
   - Window overlap-add for artifact-free output

### N-Way Extension (AudioTransportN)

Currently implements simplified N-way morphing:
- **Magnitude**: Weighted geometric mean of all sources
- **Phase**: Weighted circular mean considering both barycentric weights and local magnitudes
- **Transport**: Simplified nearest-neighbor matching (full N-way optimal transport TODO)

## Usage in AudioPlayback

### Morphing Modes:
1. **Standard Morph Mode** (`mMorphMode`):
   - Uses FluCoMa's pairwise AudioTransport
   - Activated during crossfades between audio sources
   - Smooth spectral interpolation instead of simple crossfade

2. **N-Way Morph Mode** (`mNWayMorphMode`):
   - Uses AudioTransportN for multiple source blending
   - Supports 2-8 simultaneous sources
   - Barycentric coordinate system for intuitive weight control

### Key Parameters:
- `mMorphSTFTSize`: FFT size for morphing (default 1024)
- `mMorphTransitionDuration`: Duration of morph transition in seconds
- `mMorphTargets`: Vector of (fileIndex, sampleIndex) pairs
- `mMorphWeights`: Barycentric weights for each target

## Implementation Status

### Completed:
- ✅ Basic pairwise morphing integration
- ✅ N-way barycentric coordinate system design
- ✅ AudioTransportN class with basic N-way morphing
- ✅ Thread-safe target/weight management
- ✅ Magnitude interpolation (geometric mean)
- ✅ Phase interpolation (circular mean)

### TODO:
- ⚠️ Full N-way optimal transport algorithm (currently using simplified approach)
- ⚠️ Complete ProcessNWayMorphFrame implementation in AudioPlayback
- ⚠️ SIMD optimizations for performance
- ⚠️ Von Mises distribution for improved phase interpolation
- ⚠️ UI controls for N-way morphing in ExplorerMenu

## References

- Henderson, T., & Solomon, J. (2019). Audio Transport: A Generalized Portamento via Optimal Transport. Proceedings of DAFX 2019.
- FluCoMa: Fluid Corpus Manipulation toolkit (https://www.flucoma.org/)
- Original implementation: https://github.com/sportdeath/audio_transport/