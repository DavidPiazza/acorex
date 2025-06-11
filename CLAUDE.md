# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ACorEx (3D Audio Corpus Explorer) is an openFrameworks application for exploring audio corpora using 3D visualization and time series analysis. It provides interfaces for analyzing, browsing, and listening to audio files, similar to AudioStellar and CataRT.

## Build Commands

### Initial Setup
```bash
# Clone and build all dependencies (including openFrameworks)
./build-deps.sh

# Force recompile dependencies
./build-deps.sh -c

# Force re-download and recompile
./build-deps.sh -d
```

### Building the Application
- **macOS**: Open `acorex.xcodeproj` in Xcode and build using the Debug or Release scheme
- **Windows**: Open `acorex.vcxproj` in Visual Studio and build
- **Linux**: Use `make` in the project directory

### Running Tests
No automated test framework is currently configured. Testing is manual through the application UI.

## Architecture

### Core Components

1. **Main Application** (`src/ofApp.h/cpp`)
   - Entry point and main controller
   - Manages two primary menus: AnalyserMenu and ExplorerMenu
   - Handles audio stream setup and window management

2. **Analyser System** (`src/Analyser/`)
   - `Controller`: Manages audio analysis workflow
   - `GenAnalysis`: Generic analysis framework
   - `UMAP`: Dimensionality reduction using UMAP algorithm
   - Processes audio files to extract features (pitch, loudness, shape, MFCC)

3. **Explorer System** (`src/Explorer/`)
   - `RawView`: 3D visualization of analyzed audio data points
   - `LiveView`: Real-time visualization during playback
   - `AudioPlayback`: Handles audio playback and morphing between samples
   - `PointPicker`: User interaction for selecting points in 3D space

4. **Utils** (`src/Utils/`)
   - `Data.h`: Core data structures (DataSet, AudioData, TimeData, StatsData)
   - `AudioFileLoader`: Loads various audio file formats
   - `DatasetConversion`: Converts between data formats
   - `JSON`: Handles corpus file I/O in JSON format

### Key Data Flow

1. **Analysis Phase**: Audio files → Feature extraction → Statistical reduction → JSON corpus file
2. **Exploration Phase**: Load corpus → 3D visualization → User interaction → Audio playback

### Important Configurations

- Uses C++17 (configured in `config.make`)
- Minimum macOS version: 10.15
- openFrameworks addons: ofxGui, ofxDropdown, ofxAudioFile
- External dependencies: FluCoMa, Eigen, HISSTools, foonathan/memory, Spectra, nlohmann/json

### Audio Features

- **Morph Mode**: Smooth interpolation between audio samples using STFT
- **Crossover Jump**: Probabilistic jumping between similar points in the corpus
- **Multi-dimensional navigation**: X, Y, Z axes plus color mapping to any analyzed dimension

## Development Notes

- The project structure follows openFrameworks conventions with source in `src/`
- Dependencies are built into `deps/` and `libs/` directories
- The build script handles all dependency management automatically
- Audio corpus files are stored as JSON with analyzed feature data
- The application supports both time-series and statistical views of audio data