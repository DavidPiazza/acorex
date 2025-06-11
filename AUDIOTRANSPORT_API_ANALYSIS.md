# fluid::algorithm::AudioTransport API Analysis

## Overview
AudioTransport is a spectral morphing algorithm based on optimal transport theory that creates smooth interpolations between two audio signals.

## Key Components

### Constructor
```cpp
AudioTransport(index maxFFTSize, Allocator& alloc)
```
- Creates instance with maximum FFT size
- Requires allocator for memory management

### Initialization
```cpp
void init(index windowSize, index fftSize, index hopSize)
```
- Must be called before processing
- Sets up STFT parameters
- Window size typically equals FFT size
- Hop size typically windowSize/2 or windowSize/4

### Main Processing Method
```cpp
void processFrame(RealVectorView in1, RealVectorView in2, double weight,
                  RealMatrixView out, Allocator& alloc)
```
- **in1**: First audio frame (source)
- **in2**: Second audio frame (target)
- **weight**: Interpolation weight (0.0 = all in1, 1.0 = all in2)
- **out**: Output matrix with 2 rows:
  - Row 0: Morphed audio output
  - Row 1: Window normalization values
- **alloc**: Memory allocator

### Usage Pattern
1. Create instance with max FFT size
2. Call `init()` with desired parameters
3. For each audio block:
   - Prepare input frames of size `windowSize`
   - Call `processFrame()` with interpolation weight
   - Apply window normalization to output

## Integration Strategy for AudioPlayback

### 1. Add Members to AudioPlayback
```cpp
// In AudioPlayback.h
private:
    // Morphing system
    bool mMorphMode = false;
    std::unique_ptr<fluid::algorithm::AudioTransport> mAudioTransport;
    fluid::RealMatrix mMorphOutputBuffer;
    
    // Morphing parameters
    std::atomic<int> mSTFTSize = 1024;
    std::atomic<float> mMorphDuration = 1.0f;
```

### 2. Initialize in Constructor
```cpp
// In AudioPlayback constructor
mAudioTransport = std::make_unique<fluid::algorithm::AudioTransport>(
    4096, // max FFT size
    allocator
);
mMorphOutputBuffer.resize(2, 4096);
```

### 3. Modify Crossfade Logic
Replace crossfade processing with morphing when `mMorphMode` is true:

```cpp
if (mMorphMode && mAudioTransport->initialized()) {
    // Use AudioTransport for morphing
    float weight = calculateMorphWeight(crossfadeProgress);
    
    // Get frames from both sources
    fluid::RealVector frame1(mSTFTSize);
    fluid::RealVector frame2(mSTFTSize);
    // Fill frames from audio data...
    
    // Process morph
    mAudioTransport->processFrame(
        frame1, frame2, weight,
        mMorphOutputBuffer, allocator
    );
    
    // Copy morphed output to playhead buffer
    // Apply window normalization
} else {
    // Use existing cosine/sine crossfade
}
```

### 4. Parameter Updates
- Update STFT parameters when changed (outside audio thread)
- Recalculate morph duration in samples
- Handle parameter validation

### 5. Thread Safety
- Initialize/reinitialize AudioTransport outside audio thread
- Use atomic variables for parameters
- Ensure allocator is thread-safe or use per-thread allocators

## Benefits of Integration
1. Preserves all existing AudioPlayback functionality
2. Morphing becomes an optional enhancement
3. Can switch between crossfade and morph modes
4. Leverages existing buffer management
5. Works with multiple playheads naturally