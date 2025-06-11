# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ACorEx (3D Audio Corpus Explorer) is an OpenFrameworks-based C++ application for analyzing and visualizing large audio collections. It provides 3D navigation through audio corpora with time series analysis capabilities.

## Current Development Goal: Audio Morphing Feature

We are implementing an audio morphing feature that allows smooth sonic transitions between audio files in LiveView mode, using the fluid::algorithm::AudioTransport algorithm for spectral interpolation.

### Project Status
- Task tracking via TaskMaster in `.taskmaster/tasks/tasks.json`
- 15 main tasks with 75 subtasks defined
- Starting with Task 1: Setup ofSoundStream Audio Output System

### Key Implementation Areas
1. **Audio Infrastructure**: Replace ofSoundPlayer with ofSoundStream for real-time processing
2. **UI Controls**: Add morph mode toggle and parameters to ExplorerMenu
3. **Dual Processing Modes**: Support both legacy and morph mode audio playback
4. **Buffer Management**: Use fluid::RealVector for audio data storage
5. **Spectral Morphing**: Integrate fluid::algorithm::AudioTransport for STFT-based morphing

## Build Commands

### Initial Setup
```bash
# First time setup - downloads OpenFrameworks and all dependencies
./build-deps.sh

# Force re-download dependencies
./build-deps.sh -d

# Force recompile foonathan/memory
./build-deps.sh -c
```

### Building the Project
```bash
# Build release version (optimized)
make Release

# Build debug version
make Debug

# Clean build artifacts
make clean

# Clean specific target
make CleanDebug
make CleanRelease
```

### Platform-Specific Building
- **macOS**: Open `acorex.xcodeproj` in Xcode and build
- **Windows**: Open `acorex.sln` in Visual Studio and build
- **Linux**: Use make commands (experimental support)

## Architecture Overview

The application has two main modes:
1. **Analyser Mode** - Processes audio files to extract features
2. **Explorer Mode** - Visualizes and navigates the analyzed data in 3D space

### Core Modules

**Analyser (`src/Analyser/`)**
- `Controller`: Orchestrates the analysis workflow
- `GenAnalysis`: Performs audio feature extraction (pitch, loudness, shape, MFCCs)
- `UMAP`: Dimensionality reduction for visualization

**Explorer (`src/Explorer/`)**
- `RawView`: Direct visualization of audio features
- `LiveView`: Interactive 3D navigation with audio playback
- `PointPicker`: Handle point selection in 3D space

**Utils (`src/Utils/`)**
- `Data.h`: Core data structures (TimeData, StatsData, AnalysisSettings)
- `JSON`: Serialization/deserialization of analysis data
- `DatasetConversion`: Format conversions between different data representations

### Data Flow
1. Audio files → Feature extraction (GenAnalysis)
2. Features → Statistical summaries (mean, std dev, skewness, kurtosis, quartiles)
3. High-dimensional features → UMAP reduction to 2D/3D
4. Reduced data → 3D visualization where each point is an audio file
5. User interaction → Audio playback and navigation

## Key Dependencies

- **OpenFrameworks 0.12.0**: Creative coding framework
- **ofxAudioFile**: Audio file I/O (custom fork)
- **ofxDropdown**: UI dropdown menus (custom fork)
- **flucoma-core**: Audio analysis algorithms (acorex branch)
- **Eigen 3.4.0**: Linear algebra
- **nlohmann/json**: JSON handling

## Development Notes

- C++17 standard is required
- macOS minimum version: 10.15
- No project-level testing framework (tests exist only in dependencies)
- No project-level linting configuration
- UI is built using OpenFrameworks' ofxGui addon
- Audio analysis uses modified FluCoMa algorithms

## Audio Morphing Implementation Details

### Technical Approach
1. **ofSoundStream Integration**
   - Add `ofSoundStream soundStream` to ofApp
   - Implement `audioOut(float* output, int bufferSize, int nChannels)`
   - Route audio from LiveView via `getAudioBlock()`

2. **Morph Mode UI (ExplorerMenu)**
   - `ofxToggle mMorphModeToggle` - Enable/disable morphing
   - `ofxFloatSlider mTransitionDurationSlider` - Control morph duration (0.1-5.0s)
   - `ofxIntSlider mSTFTSizeSlider` - STFT window size (512-4096)
   - `ofxFloatSlider mInterpolationSlider` - Manual morph control
   - `ofxToggle mLoopToggle` - Loop morphing

3. **LiveView Audio Processing**
   - Dual mode support: Legacy (ofSoundPlayer) vs Morph (AudioTransport)
   - Buffer management with `fluid::RealVector mCurrentAudioBuffer`, `mTargetAudioBuffer`
   - Interpolation weight management (0.0 to 1.0 over transition duration)
   - Single sound passthrough when only one file selected

4. **fluid::algorithm::AudioTransport**
   - STFT-based spectral morphing between audio streams
   - Parameters: Window size, hop size (window/2), FFT size (=window)
   - `processFrame()` called per audio block with interpolation weight

### Current Implementation Focus
Working on Task 1: Setting up the ofSoundStream infrastructure as the foundation for all audio processing.

## CRITICAL UPDATE: Architecture Mismatch Discovery

We discovered that the original v1.0.0 release has a completely different audio architecture:

### Reference Code (v1.0.0) Architecture:
- **AudioPlayback class**: Sophisticated polyphonic audio engine with:
  - Multiple simultaneous playheads
  - Probabilistic jump system
  - Crossfading between segments  
  - Device/buffer selection
- **AudioFileLoader class**: Dedicated audio loading with resampling
- **Additional UI controls**: Jump chances, crossfade settings, device selection, etc.

### Our Current Implementation:
- Built on a fork that removed the AudioPlayback system
- Replaced everything with ofSoundPlayer (legacy)
- Added our own ofSoundStream + AudioTransport morphing on top

### The Problem:
Our morphing implementation bypassed the existing sophisticated audio system. We need to integrate morphing INTO the AudioPlayback system, not replace it.

### New Plan:
1. **STOP current implementation** - It's built on wrong assumptions
2. **Copy reference code files** - Restore AudioPlayback and AudioFileLoader 
3. **Integrate morphing into AudioPlayback** - Add spectral morphing as a feature of the existing system
4. **Preserve existing features** - Keep probabilistic jumps, crossfading, etc.

### Integration Strategy:
The AudioPlayback system uses a callback `audioOut()` that we can enhance:
- When morph mode is OFF: Use existing crossfade/jump behavior
- When morph mode is ON: Apply AudioTransport spectral morphing between playheads
- Morph parameters control the spectral interpolation instead of position jumps