# AudioPlayback System Architecture Analysis

## Overview
The AudioPlayback class is a sophisticated polyphonic audio engine with the following key features:

### Core Architecture
1. **Multiple Playheads**: Supports multiple simultaneous audio playheads, each playing from different positions in the corpus
2. **Probabilistic Jumps**: Playheads can probabilistically jump to nearby points in the audio corpus based on spatial distance
3. **Crossfading**: Smooth transitions between jump points using cosine/sine crossfade curves
4. **Thread-Safe Design**: Uses mutexes and atomic variables for safe multi-threaded operation

### Key Components

#### AudioPlayhead Structure
```cpp
struct AudioPlayhead {
    size_t playheadID;
    size_t fileIndex;       // Current file being played
    size_t sampleIndex;     // Current position in file
    
    // Crossfading state
    bool crossfading;
    size_t jumpFileIndex;   // Target file for crossfade
    size_t jumpSampleIndex; // Target position
    size_t crossfadeCurrentSample;
    size_t crossfadeSampleLength;
    
    std::queue<size_t> triggerSamplePoints;  // Points where jumps can occur
};
```

#### Audio Processing Pipeline (audioOut method)
1. **Playhead Management**: Add new playheads and remove finished ones
2. **Per-Playhead Processing**:
   - Handle crossfading if active
   - Check trigger points for potential jumps
   - Fill audio buffer segments
   - Apply probabilistic jump logic
   - Handle end-of-file (loop or kill)
3. **Buffer Mixing**: Mix all playhead outputs into final stereo buffer

#### Jump System
- **Trigger Points**: Calculated based on FFT hop size (analysis frame boundaries)
- **Jump Probability**: Controlled by `mCrossoverJumpChanceX1000`
- **Spatial Constraints**: 
  - Maximum jump distance in 3D space
  - Maximum number of jump targets to consider
  - Option to allow/disallow jumps within same file
  - Minimum time difference for same-file jumps

### Integration Points for Morphing

1. **Replace Crossfade with Morph Option**:
   - Current: Uses cosine/sine crossfade between jump points
   - Enhancement: Add option to use AudioTransport morphing instead

2. **Add Morph Mode Flag**:
   - Add `bool mMorphMode` to control behavior
   - When false: Use existing crossfade
   - When true: Use AudioTransport morphing

3. **AudioTransport Integration**:
   - Add `fluid::algorithm::AudioTransport mAudioTransport` member
   - Initialize with appropriate STFT parameters
   - Process frames during crossfade/morph sections

4. **Parameter Management**:
   - Morph duration (replace/supplement crossfade length)
   - STFT window size
   - Interpolation control

5. **Buffer Management**:
   - Current system already handles dual buffers during crossfade
   - Can reuse same infrastructure for morphing

### Thread Safety Considerations
- Audio callback runs on separate thread
- Must maintain lock-free audio processing where possible
- Parameter updates should use atomic variables
- AudioTransport initialization should happen outside audio thread

### Next Steps
1. Add AudioTransport member and initialization
2. Create morph mode flag and parameters
3. Modify crossfade logic to optionally use morphing
4. Connect UI controls to parameters
5. Test integration with existing features