# AudioTransport Morphing Integration Design

## Architecture Overview

The morphing feature will be integrated INTO the existing AudioPlayback system as an optional enhancement that works alongside the current probabilistic jump and crossfade features.

## Design Principles

1. **Preserve Existing Functionality**: All current features must continue to work
2. **Optional Enhancement**: Morphing is an alternative to crossfading, not a replacement
3. **Minimal Disruption**: Changes should be localized to specific areas
4. **Thread Safety**: Maintain existing thread-safe design patterns

## Integration Points

### 1. AudioPlayback Class Modifications

#### New Member Variables
```cpp
// In AudioPlayback.h private section:

// Morphing system
std::atomic<bool> mMorphMode{false};
std::unique_ptr<fluid::algorithm::AudioTransport> mAudioTransport;
fluid::RealMatrix mMorphOutputBuffer;
Allocator mAllocator;

// Morphing parameters
std::atomic<int> mMorphSTFTSize{1024};
std::atomic<float> mMorphTransitionDuration{1.0f};
std::atomic<bool> mMorphParametersChanged{false};

// Thread-safe parameter staging
std::mutex mMorphParamMutex;
int mStagedSTFTSize = 1024;
```

#### New Methods
```cpp
// Public methods for UI control
void SetMorphMode(bool enabled);
void SetMorphSTFTSize(int size);
void SetMorphTransitionDuration(float duration);

// Private methods
void UpdateMorphParameters(); // Called at start of audioOut
void ProcessMorphFrame(ofSoundBuffer* outBuffer, size_t* outBufferPosition, 
                      Utils::AudioPlayhead* playhead, size_t crossfadeLength);
```

### 2. Modified Audio Processing Flow

#### In audioOut() method:
```cpp
// At the beginning of audioOut
if (mMorphParametersChanged) {
    UpdateMorphParameters(); // Safe parameter update
}

// In the crossfade processing section:
if (mPlayheads[playheadIndex].crossfading) {
    if (mMorphMode && mAudioTransport && mAudioTransport->initialized()) {
        ProcessMorphFrame(&playheadBuffer, &playheadBufferPosition, 
                         &mPlayheads[playheadIndex], crossfadeSamplesLeft);
    } else {
        // Existing cosine/sine crossfade code
    }
}
```

### 3. Morph Processing Implementation

```cpp
void AudioPlayback::ProcessMorphFrame(ofSoundBuffer* outBuffer, 
                                     size_t* outBufferPosition,
                                     Utils::AudioPlayhead* playhead, 
                                     size_t crossfadeLength) {
    // Calculate interpolation weight based on crossfade progress
    float weight = (float)playhead->crossfadeCurrentSample / 
                   (float)playhead->crossfadeSampleLength;
    
    // Prepare input frames
    size_t frameSize = mMorphSTFTSize;
    fluid::RealVector frame1(frameSize);
    fluid::RealVector frame2(frameSize);
    
    // Fill frames from audio data with bounds checking
    for (size_t i = 0; i < frameSize; i++) {
        size_t idx1 = playhead->sampleIndex + i;
        size_t idx2 = playhead->jumpSampleIndex + i;
        
        frame1[i] = (idx1 < audioSize1) ? 
            mRawView->GetAudioData()->raw[playhead->fileIndex].getSample(idx1, 0) : 0.0;
        frame2[i] = (idx2 < audioSize2) ? 
            mRawView->GetAudioData()->raw[playhead->jumpFileIndex].getSample(idx2, 0) : 0.0;
    }
    
    // Process morph
    mAudioTransport->processFrame(frame1, frame2, weight, 
                                 mMorphOutputBuffer, mAllocator);
    
    // Copy morphed output with window normalization
    auto morphed = mMorphOutputBuffer.row(0);
    auto window = mMorphOutputBuffer.row(1);
    
    for (size_t i = 0; i < crossfadeLength && i < frameSize; i++) {
        float sample = morphed[i] / (window[i] > 0 ? window[i] : 1.0);
        outBuffer->getSample(*outBufferPosition + i, 0) = sample;
    }
    
    // Update positions
    playhead->crossfadeCurrentSample += crossfadeLength;
    playhead->sampleIndex += crossfadeLength;
    playhead->jumpSampleIndex += crossfadeLength;
    *outBufferPosition += crossfadeLength;
}
```

### 4. UI Integration

#### ExplorerMenu Connections
```cpp
// In ExplorerMenu::Initialise()
mMorphModeToggle.addListener(this, &ExplorerMenu::onMorphModeChanged);
mSTFTSizeSlider.addListener(this, &ExplorerMenu::onSTFTSizeChanged);
mTransitionDurationSlider.addListener(this, &ExplorerMenu::onTransitionDurationChanged);

// Callback implementations
void ExplorerMenu::onMorphModeChanged(bool& enabled) {
    mLiveView.GetAudioPlayback()->SetMorphMode(enabled);
}

void ExplorerMenu::onSTFTSizeChanged(int& size) {
    mLiveView.GetAudioPlayback()->SetMorphSTFTSize(size);
}

void ExplorerMenu::onTransitionDurationChanged(float& duration) {
    mLiveView.GetAudioPlayback()->SetMorphTransitionDuration(duration);
}
```

### 5. Parameter Management

#### Thread-Safe Updates
```cpp
void AudioPlayback::SetMorphSTFTSize(int size) {
    std::lock_guard<std::mutex> lock(mMorphParamMutex);
    mStagedSTFTSize = size;
    mMorphParametersChanged = true;
}

void AudioPlayback::UpdateMorphParameters() {
    std::lock_guard<std::mutex> lock(mMorphParamMutex);
    if (mStagedSTFTSize != mMorphSTFTSize) {
        mMorphSTFTSize = mStagedSTFTSize;
        // Reinitialize AudioTransport with new size
        mAudioTransport->init(mMorphSTFTSize, mMorphSTFTSize, mMorphSTFTSize/2);
    }
    mMorphParametersChanged = false;
}
```

### 6. Memory Management

- Use OpenFrameworks allocator or create custom allocator
- Pre-allocate buffers in constructor
- Ensure all allocations happen outside audio thread

### 7. Transition Duration Handling

The morph transition duration will map to the existing crossfade length:
```cpp
// When morph mode is enabled and a jump is triggered:
if (mMorphMode) {
    // Convert duration in seconds to samples
    int morphSamples = mMorphTransitionDuration * mSoundStream.getSampleRate();
    mPlayheads[playheadIndex].crossfadeSampleLength = morphSamples;
}
```

## Benefits of This Design

1. **Minimal Code Changes**: Most changes are additions rather than modifications
2. **Backward Compatible**: Existing functionality remains unchanged
3. **Performance**: Only active when morph mode is enabled
4. **Flexible**: Can easily switch between modes
5. **Thread Safe**: Follows existing patterns for parameter updates

## Testing Strategy

1. **Unit Tests**: Test AudioTransport integration in isolation
2. **Integration Tests**: Verify morphing works with multiple playheads
3. **Regression Tests**: Ensure crossfade mode still works correctly
4. **Performance Tests**: Measure CPU usage with morphing enabled
5. **UI Tests**: Verify parameter changes are applied correctly