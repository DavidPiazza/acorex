# CLAUDE.md - Utils/

This file provides guidance to Claude Code (claude.ai/code) when working with code in the Utils directory.

## Module Overview

The Utils module provides core data structures, file I/O operations, and shared definitions used throughout the ACorEx application.

## Core Components

### Data.h
Central data structure definitions for the entire application.

**Primary Structures:**
- **DataSet**: Main container holding all corpus data
  - `dimensionNames`: Feature names (e.g., "Pitch Mean", "MFCC_1 Std Dev")
  - `fileList`: Paths to audio files in the corpus
  - `audio`: Audio buffers (loaded on demand)
  - `time`: Time-series analysis data
  - `stats`: Statistical summaries
  - `analysisSettings`: Configuration used during analysis

- **AudioData**: Audio buffer storage
  - `loaded`: Flags indicating which files are loaded
  - `raw`: ofSoundBuffer objects containing audio samples

- **TimeData**: Frame-by-frame analysis results
  - 3D array: [file][timepoint][dimension]
  - First dimension is always time in seconds

- **StatsData**: Statistical summaries
  - `raw`: [file][dimension][statistic] for unreduced data
  - `reduced`: [file][dimension] for UMAP-reduced data

- **Playhead Structures**: State tracking for audio playback
  - `AudioPlayhead`: Audio engine state
  - `VisualPlayhead`: Visualization state

### JSON.h/cpp
Serialization and deserialization of corpus data.

**Key Functions:**
- `readDataSet()`: Loads corpus from JSON file
- `writeDataSet()`: Saves corpus to JSON file
- `packSettings()`: Serializes analysis configuration
- `unpackSettings()`: Deserializes analysis configuration

**JSON Structure:**
```json
{
  "analysisSettings": { ... },
  "dimensionNames": [ ... ],
  "fileList": [ ... ],
  "stats": {
    "raw": [ ... ],
    "reduced": [ ... ]
  },
  "time": {
    "raw": [ ... ]
  }
}
```

### AudioFileLoader.h/cpp
Audio file reading and preprocessing.

**Features:**
- Supports multiple formats via ofxAudioFile
- Automatic mono conversion
- Resampling to target sample rate
- Returns FluCoMa-compatible RealVector format
- Memory efficient loading

**Key Methods:**
- `readFile()`: Main loading function
- `convertToMono()`: Stereo to mono conversion
- `resample()`: Sample rate conversion

### DatasetConversion.h/cpp
Utilities for converting between different data representations.

**Conversions:**
- Time-series to statistics
- Statistics to flat arrays for UMAP
- Dimension index mapping
- Corpus format migrations

### DimensionBounds.h
Helper for calculating min/max bounds across dimensions.

**Usage:**
- Normalizing data for visualization
- Color mapping ranges
- Axis scaling in 3D view

### InterfaceDefs.h
Shared UI constants and styling definitions.

**Includes:**
- Color schemes (`Colors` struct)
- Layout constants (`MenuLayout` struct)
- UI element sizes and spacing
- Font settings

## Design Patterns

### Data Flow
1. **Loading**: JSON file → JSON::readDataSet → DataSet populated
2. **Access**: Components request data via typed accessors
3. **Saving**: DataSet → JSON::writeDataSet → JSON file

### Memory Management
- Audio data loaded on-demand to conserve memory
- Shared pointers for DataSet to enable safe sharing
- RAII patterns for resource management

### Thread Safety
- DataSet is treated as immutable after loading
- No direct modification during exploration
- New corpus creation generates new DataSet instance

## Integration Notes

- All modules depend on Utils for data structures
- JSON format enables corpus portability
- Audio loading uses both ofxAudioFile and FluCoMa
- Consistent error handling with descriptive messages