# CLAUDE.md - Explorer/

This file provides guidance to Claude Code (claude.ai/code) when working with code in the Explorer directory.

## Module Overview

The Explorer module enables interactive exploration and playback of analyzed audio corpora through 3D visualization and sophisticated audio morphing/crossfading capabilities.

## Core Components

### AudioPlayback.h/cpp
The main audio engine that handles real-time playback and morphing between corpus points.

**Key Features:**
- Multiple simultaneous playheads with independent control
- Two playback modes:
  - **Morph Mode**: Smooth spectral interpolation between samples using FluCoMa AudioTransport
  - **Crossover Mode**: Probabilistic jumping with crossfades between similar points
- Thread-safe audio generation for real-time performance
- Integration with PointPicker for user-driven selection
- Visual playhead tracking for UI feedback

**Important Methods:**
- `getBuffer()`: Main audio callback, generates audio samples
- `setPlayheadTargets()`: Updates morph destinations
- `checkForJumps()`: Implements crossover jump logic
- `updateVisualPlayheads()`: Syncs visual representation with audio state

### RawView.h/cpp
Data access layer for corpus exploration.

**Responsibilities:**
- Loading corpus data from JSON files
- Providing typed access to time/stats/audio data
- Checking corpus properties (analysis types, reduction status)
- Managing the shared DataSet instance

**Key Methods:**
- `loadCorpus()`: Loads and validates corpus files
- `getTimeDataRef()`: Access time-series analysis data
- `getStatsDataRef()`: Access statistical analysis data
- `dimensionIsTime()`: Check if corpus includes time analysis

### LiveView.h/cpp
Real-time 3D visualization of the corpus and playheads.

**Features:**
- Interactive camera controls
- Point cloud rendering with dimension-based coloring
- Playhead visualization and trails
- Picking support for point selection

### PointPicker.h/cpp
Handles user interaction for selecting points in 3D space.

**Functionality:**
- Ray-casting for 3D point selection
- Nearest-neighbor search
- Integration with camera system
- Selection feedback

### SpaceDefs.h
Common definitions for spatial calculations and constants used across Explorer components.

## Data Flow

1. **Corpus Loading**: RawView loads JSON → validates → provides data access
2. **User Interaction**: Mouse/keyboard → PointPicker → selected points
3. **Audio Generation**: Selected points → AudioPlayback → morphed/crossfaded audio
4. **Visualization**: Corpus data + playhead positions → LiveView → 3D rendering

## Threading Model

- **Main Thread**: UI, visualization, user input
- **Audio Thread**: AudioPlayback::getBuffer() runs on audio callback thread
- **Synchronization**: Lock-free audio generation with atomic operations for parameter updates

## Integration Points

- Receives corpus data from Analyser module via JSON files
- Uses Utils module data structures (DataSet, AudioData, etc.)
- Controlled by ExplorerMenu for UI integration
- Audio output routed through ofApp audio callbacks