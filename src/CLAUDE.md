# CLAUDE.md - src/

This file provides guidance to Claude Code (claude.ai/code) when working with code in the src directory.

## Directory Overview

The `src/` directory contains the main application code for ACorEx, organized into three primary modules plus the application framework files.

## Main Application Files

- **main.cpp**: Entry point, sets up OpenGL window and launches the app
- **ofApp.h/cpp**: Main application controller that:
  - Manages the two primary UI menus (AnalyserMenu and ExplorerMenu)
  - Handles audio stream setup and routing
  - Manages window events and UI layout
  - Coordinates switching between analysis and exploration modes

## Module Organization

### 1. Analyser/
Audio corpus creation and analysis module. Processes audio files to extract features and perform dimensionality reduction.

### 2. Explorer/
Interactive exploration module. Provides 3D visualization and playback of analyzed audio corpora.

### 3. Utils/
Shared utilities, data structures, and file I/O operations used across the application.

## Key Interactions

- **AnalyserMenu.h/cpp**: UI wrapper for the Analyser module
  - Manages analysis settings and file selection
  - Triggers corpus creation and reduction operations
  - Handles progress feedback during processing

- **ExplorerMenu.h/cpp**: UI wrapper for the Explorer module  
  - Manages corpus loading and dimension selection
  - Controls playback modes (morph vs crossover)
  - Handles audio device configuration
  - Coordinates between RawView (data) and LiveView (visualization)

## Audio Flow

1. Audio output is generated in `Explorer/AudioPlayback`
2. ExplorerMenu::getAudioBlock() retrieves audio from AudioPlayback
3. ofApp::audioOut() sends this to the audio device via ofSoundStream

## State Management

- The application maintains separate states for analysis and exploration
- Only one mode can be active at a time (controlled by toggle buttons)
- Audio streaming is only active during exploration mode