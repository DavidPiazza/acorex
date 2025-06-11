# CLAUDE.md - Analyser/

This file provides guidance to Claude Code (claude.ai/code) when working with code in the Analyser directory.

## Module Overview

The Analyser module processes audio files to create searchable corpora by extracting features and optionally reducing dimensionality for visualization.

## Core Components

### Controller.h/cpp
Main orchestrator for all analysis operations.

**Primary Functions:**
- `CreateCorpus()`: Analyzes a directory of audio files to create a new corpus
- `ReduceCorpus()`: Applies UMAP dimensionality reduction to high-dimensional features
- `InsertIntoCorpus()`: Adds new audio files to an existing corpus
- `MergeDatasets()`: Combines multiple corpora while preserving analysis settings

**Key Responsibilities:**
- Directory traversal and audio file discovery
- Progress tracking and cancellation support
- Coordination between GenAnalysis (feature extraction) and UMAP (reduction)
- Dimension naming based on analysis settings
- Thread management for parallel processing

### GenAnalysis.h/cpp
Generic audio analysis framework for feature extraction.

**Supported Features:**
- **Pitch**: Fundamental frequency estimation
- **Loudness**: Perceptual loudness analysis
- **Shape**: Spectral shape descriptors
- **MFCC**: Mel-frequency cepstral coefficients
- **Time-series**: Frame-by-frame analysis
- **Statistics**: Mean, std deviation, skewness, kurtosis, percentiles

**Analysis Pipeline:**
1. Load audio file → mono conversion → resampling
2. Apply windowing (FFT size, hop size)
3. Extract features per frame
4. Calculate statistics across frames
5. Store in DataSet structure

### UMAP.h/cpp
Dimensionality reduction using Uniform Manifold Approximation and Projection.

**Functionality:**
- Reduces high-dimensional feature space to 2D/3D for visualization
- Preserves local neighborhood structure
- Configurable parameters (dimensions, iterations)
- Progress callback support
- Integration with statistical data from GenAnalysis

## Analysis Workflow

1. **Corpus Creation**:
   ```
   Audio files → GenAnalysis → Feature extraction → Statistics → JSON output
   ```

2. **Dimensionality Reduction**:
   ```
   High-dim corpus → UMAP → Low-dim representation → Updated JSON
   ```

3. **Corpus Expansion**:
   ```
   Existing corpus + New files → Merged analysis → Consistent dimensions
   ```

## Configuration Parameters

**Analysis Settings** (from Utils/Data.h):
- `windowFFTSize`: FFT window size (default: 1024)
- `hopFraction`: Hop size as fraction of window (default: 2)
- `nBands`: Number of frequency bands for shape analysis
- `nCoefs`: Number of MFCC coefficients
- `minFreq`/`maxFreq`: Frequency range for analysis
- `sampleRate`: Target sample rate for resampling

## Progress Tracking

All operations support progress callbacks:
- File-level progress during analysis
- Iteration progress during UMAP reduction
- Cancellation via boolean flag reference

## Error Handling

- Validates audio file formats before processing
- Checks corpus compatibility for merging operations
- Handles missing files gracefully
- Reports specific errors through console output

## Integration Notes

- Output format is JSON (handled by Utils/JSON)
- Uses FluCoMa for audio analysis algorithms
- Parallel processing with configurable thread count
- Memory efficient streaming for large corpora