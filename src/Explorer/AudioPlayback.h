/*
The MIT License (MIT)

Copyright (c) 2024 Elowyn Fearne

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include "./RawView.h"
#include "./PointPicker.h"
#include "Utils/Data.h"
#include <ofSoundBuffer.h>
#include <ofSoundStream.h>
#include <ofMesh.h>
#include <vector>
#include <mutex>
#include <atomic>
#include <memory>
#include <algorithms/public/AudioTransport.hpp>
#include "AudioTransportN.hpp"
#include "AudioMorphEngine.h"
#include <data/TensorTypes.hpp>
#include <data/FluidMemory.hpp>

namespace Acorex {
namespace Explorer {

class AudioPlayback {
public:
	// Playback mode enumeration
	enum class PlaybackMode {
		DISCRETE,    // Traditional jump-based playback with crossfades
		CONTINUOUS   // Continuous morphing through AudioMorphEngine
	};

	AudioPlayback ( ) { }
	~AudioPlayback ( ) { }

	void Initialise ( );
	void RestartAudio ( size_t sampleRate, size_t bufferSize, ofSoundDevice outDevice );

	void audioOut ( ofSoundBuffer& outBuffer );

	void SetRawView ( std::shared_ptr<RawView>& rawPointer ) { mRawView = rawPointer; }

	bool CreatePlayhead ( size_t fileIndex, size_t sampleIndex );
	bool KillPlayhead ( size_t playheadID );
	std::vector<Utils::VisualPlayhead> GetPlayheadInfo ( );
	void SetFlagReset ( );
	void WaitForResetConfirm ( );

	void SetTimeCorpus ( const std::vector<ofMesh>& timeCorpus );

	void SetPointPicker ( std::shared_ptr<PointPicker>& pointPicker ) { mPointPicker = pointPicker; }

	void SetLoopPlayheads ( bool loop ) { mLoopPlayheads = loop; }
	void SetJumpSameFileAllowed ( bool allowed ) { mJumpSameFileAllowed = allowed; }
	void SetJumpSameFileMinTimeDiff ( int timeDiff ) { mJumpSameFileMinTimeDiff = timeDiff; }
	void SetCrossoverJumpChance ( int jumpsInAThousand ) { mCrossoverJumpChanceX1000 = jumpsInAThousand; }
	void SetCrossfadeSampleLength ( int length ) { mCrossfadeSampleLength = length; }
	void SetMaxJumpDistanceSpace ( int distanceX1000 ) { mMaxJumpDistanceSpaceX1000 = distanceX1000; }
	void SetMaxJumpTargets ( int targets ) { mMaxJumpTargets = targets; }
	
	// Morphing controls
	void SetMorphMode ( bool enabled ) { mMorphMode = enabled; }
	void SetMorphSTFTSize ( int size );
	void SetMorphTransitionDuration ( float duration ) { mMorphTransitionDuration = duration; }
	
	// N-way morphing controls
	void SetNWayMorphMode ( bool enabled ) { mNWayMorphMode = enabled; }
	void SetMorphTargets ( const std::vector<std::pair<size_t, size_t>>& targets,
	                      const AudioTransportN::BarycentricWeights& weights );
	void ClearMorphTargets ( );
	bool SetMorphTargetsFromKNN ( const glm::vec3& position, int k = 3 );

	int GetSampleRate() const { return mSampleRate.load(); }
	
	// Playback mode controls
	void SetPlaybackMode ( PlaybackMode mode );
	PlaybackMode GetPlaybackMode ( ) const { return mPlaybackMode.load(); }
	bool IsDiscreteMode ( ) const { return mPlaybackMode.load() == PlaybackMode::DISCRETE; }
	bool IsContinuousMode ( ) const { return mPlaybackMode.load() == PlaybackMode::CONTINUOUS; }
	
	// AudioMorphEngine integration
	std::shared_ptr<AudioMorphEngine> GetAudioMorphEngine ( ) { return mAudioMorphEngine; }
	void SetKDTreeForMorphEngine ( const std::shared_ptr<fluid::algorithm::KDTree>& kdTree );
	
	// Sub-frame accurate playhead tracking
	double GetPlayheadPosition ( size_t playheadID ) const;
	std::vector<std::pair<size_t, double>> GetAllPlayheadPositions ( ) const;
	
	// Mode transition controls
	void SetModeTransitionDuration ( float seconds ) { mModeTransitionDuration = seconds; }
	float GetModeTransitionDuration ( ) const { return mModeTransitionDuration.load(); }
	bool IsModeTransitioning ( ) const { return mModeTransitioning.load(); }

private:

	void FillAudioSegment ( ofSoundBuffer* outBuffer, size_t* outBufferPosition, Utils::AudioPlayhead* playhead, bool outBufferFull );
	void CrossfadeAudioSegment ( ofSoundBuffer* outBuffer, size_t* outBufferPosition, size_t startSample_A, size_t endSample_A, size_t fileIndex_A, Utils::AudioPlayhead* playhead_B, size_t lengthSetting, bool outBufferFull );

	void CalculateTriggerPoints ( Utils::AudioPlayhead& playhead );
	bool LoadAudioFile ( size_t fileIndex );
	
	// Morphing methods
	void UpdateMorphParameters ( );
	void ProcessMorphFrame ( ofSoundBuffer* outBuffer, size_t* outBufferPosition, 
	                        Utils::AudioPlayhead* playhead, size_t crossfadeLength );
	void ProcessNWayMorphFrame ( ofSoundBuffer* outBuffer, size_t* outBufferPosition,
	                            Utils::AudioPlayhead* playhead, size_t crossfadeLength );
	
	// Mode transition processing
	void ProcessModeTransition ( ofSoundBuffer& outBuffer );
	void ProcessDiscreteMode ( ofSoundBuffer& outBuffer );

	std::vector<Utils::AudioPlayhead> mPlayheads;

	std::shared_ptr<RawView> mRawView;
	std::shared_ptr<PointPicker> mPointPicker;

	ofSoundStream mSoundStream;

	bool bStreamStarted = false;

	// settings -----------------------------------
	
	std::atomic<bool> mLoopPlayheads = false;
	std::atomic<bool> mJumpSameFileAllowed = false;
	std::atomic<int> mJumpSameFileMinTimeDiff = 2;
	std::atomic<int> mCrossoverJumpChanceX1000 = 50;
	std::atomic<int> mCrossfadeSampleLength = 256;
	std::atomic<int> mMaxJumpDistanceSpaceX1000 = 50;
	std::atomic<int> mMaxJumpTargets = 5;

	// thread safety ------------------------------

	std::atomic<int> mActivePlayheads = 0;

	std::mutex mNewPlayheadMutex;
	std::queue<Utils::AudioPlayhead> mNewPlayheads;
	std::queue<size_t> mPlayheadsToKill;
	size_t playheadCounter = 0;

	std::mutex mVisualPlayheadUpdateMutex;
	std::vector<Utils::VisualPlayhead> mVisualPlayheads;

	std::mutex mTimeCorpusMutex;
	std::vector<ofMesh> mTimeCorpus;

	std::atomic<bool> bResetFlag = false;
	
	// Morphing system ----------------------------
	
	std::atomic<bool> mMorphMode{false};
	std::unique_ptr<fluid::algorithm::AudioTransport> mAudioTransport;
	std::unique_ptr<AudioTransportN> mAudioTransportN;
	fluid::RealMatrix mMorphOutputBuffer;
	fluid::Allocator mAllocator;
	
	// N-way morphing
	std::atomic<bool> mNWayMorphMode{false};
	std::mutex mMorphTargetsMutex;
	std::vector<std::pair<size_t, size_t>> mMorphTargets; // (fileIndex, sampleIndex) pairs
	AudioTransportN::BarycentricWeights mMorphWeights;
	
	// Morphing parameters
	std::atomic<int> mMorphSTFTSize{1024};
	std::atomic<float> mMorphTransitionDuration{1.0f};
	std::atomic<bool> mMorphParametersChanged{false};
	
	// Thread-safe parameter staging
	std::mutex mMorphParamMutex;
	int mStagedSTFTSize = 1024;

	// ---------------- Audio stream replacement (single global stream) ----------------
	// The playback engine no longer owns its own ofSoundStream. We keep the current
	// sample-rate so that internal calculations remain correct.
	std::atomic<int> mSampleRate{44100};
	std::atomic<size_t> mBufferSize{512};
	
	// ---------------- Playback mode management ----------------
	std::atomic<PlaybackMode> mPlaybackMode{PlaybackMode::DISCRETE};
	std::shared_ptr<AudioMorphEngine> mAudioMorphEngine;
	
	// Configuration flags for mode behavior
	bool mDefaultToContinuousMode = false;
	
	// Mode transition state
	std::atomic<bool> mModeTransitioning{false};
	std::atomic<float> mModeTransitionProgress{0.0f};
	std::atomic<float> mModeTransitionDuration{0.5f}; // seconds
	PlaybackMode mTransitionTargetMode{PlaybackMode::DISCRETE};
	ofSoundBuffer mTransitionBuffer;
	std::mutex mTransitionMutex;
};

} // namespace Explorer
} // namespace Acorex