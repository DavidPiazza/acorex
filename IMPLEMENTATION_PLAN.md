# Audio Morphing Implementation Plan - REVISED

## Current Situation
We discovered that our current implementation was built on a fork that removed the sophisticated AudioPlayback system. We need to restore the original architecture and integrate morphing into it, not replace it.

## Step-by-Step Implementation Plan

### Phase 1: Restore Missing Components
1. **Copy AudioPlayback files**
   - Copy `Explorer/AudioPlayback.cpp/h` from reference code
   - Copy `Utils/AudioFileLoader.cpp/h` from reference code
   
2. **Update LiveView to use AudioPlayback**
   - Restore the `AudioPlayback mAudioPlayback` member
   - Remove our added audio processing code
   - Restore CreatePlayhead/KillPlayhead functionality
   
3. **Restore missing UI controls in ExplorerMenu**
   - Loop playheads toggle
   - Jump settings (same file allowed, min time diff)
   - Crossover jump chance
   - Crossfade settings
   - Max jump distance/targets
   - Buffer size dropdown
   - Output device dropdown

### Phase 2: Integrate Morphing into AudioPlayback
1. **Add AudioTransport to AudioPlayback class**
   - Add `fluid::algorithm::AudioTransport` member
   - Add morph mode flag and parameters
   
2. **Modify AudioPlayback::audioOut()**
   - When morph mode OFF: Use existing behavior
   - When morph mode ON: Apply spectral morphing between playheads
   
3. **Add Morph UI Controls**
   - Keep all existing controls
   - Add morph panel with:
     - Morph mode toggle
     - STFT size
     - Transition duration
     - Interpolation position

### Phase 3: Integration Logic
1. **Morph between playheads**
   - Use existing playhead system
   - Morph between active playheads instead of jumping
   - Transition duration controls morph speed
   
2. **Preserve existing features**
   - Probabilistic jumps still work (but morph to target)
   - Crossfading still available as alternative
   - All spatial constraints preserved

## Key Differences from Current Implementation
- We're ADDING morphing to AudioPlayback, not replacing it
- All existing features remain functional
- Morphing becomes an enhancement option, not the only mode
- User can toggle between jump/crossfade and morph behaviors

## Next Immediate Steps
1. Copy the 4 missing files from reference code
2. Fix compilation errors
3. Test that original functionality works
4. Then add morphing as an enhancement