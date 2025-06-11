#include "./AudioPlayback.h"
#include "ofLog.h"
#include <random>
#include <cmath>
#include <climits>
#include <chrono>
#include <ofxAudioFile.h>
#include "Utils/AudioFileLoader.h"

using namespace Acorex;

void Explorer::AudioPlayback::Initialise ( )
{
	ofLogNotice ( "AudioPlayback::Initialise" ) << "Called";
	srand ( time ( NULL ) );
	
	// Initialize morphing system
	try {
		// Initialize standard AudioTransport
		mAudioTransport = std::make_unique<fluid::algorithm::AudioTransport>(4096, mAllocator);
		mMorphOutputBuffer = fluid::RealMatrix(2, 4096);
		mAudioTransport->init(mMorphSTFTSize, mMorphSTFTSize, mMorphSTFTSize / 2);
		
		// Initialize N-way AudioTransport
		mAudioTransportN = std::make_unique<AudioTransportN>(4096, mAllocator);
		mAudioTransportN->initN(mMorphSTFTSize, mMorphSTFTSize, mMorphSTFTSize / 2);
		
		ofLogNotice ( "AudioPlayback" ) << "Morphing systems initialized successfully";
	} catch ( const std::exception& e ) {
		ofLogError ( "AudioPlayback" ) << "Failed to initialize AudioTransport: " << e.what();
		mAudioTransport.reset();
		mAudioTransportN.reset();
		mMorphMode = false;
		mNWayMorphMode = false;
	}

	ofSoundDevice outDevice;
	std::vector<ofSoundDevice> devices = mSoundStream.getDeviceList ( ofSoundDevice::Api::DEFAULT );
	
	// Find the default output device
	for ( auto& device : devices )
	{
		if ( device.isDefaultOutput && device.outputChannels > 0 )
		{
			outDevice = device;
			break;
		}
	}
	
	// If no default found, use first available output device
	if ( outDevice.deviceID == -1 )
	{
		for ( auto& device : devices )
		{
			if ( device.outputChannels > 0 )
			{
				outDevice = device;
				break;
			}
		}
	}
	
	RestartAudio ( 44100, 512, outDevice );
}

void Explorer::AudioPlayback::RestartAudio ( size_t sampleRate, size_t bufferSize, ofSoundDevice /*outDevice*/ )
{
	ofLogNotice ( "AudioPlayback::RestartAudio" ) << "(processor only) sampleRate=" << sampleRate << " bufferSize=" << bufferSize;
	
	// Update cached stream information for later calculations
	mSampleRate = static_cast<int>( sampleRate );
	mBufferSize = bufferSize;

	// Signal currently playing heads (if any) to reset so they don't access buffers sized for the old stream
	if ( bStreamStarted )
	{
		SetFlagReset();
		WaitForResetConfirm();
	}

	bStreamStarted = true; // Processor is considered ready once parameters are set
}

void Explorer::AudioPlayback::audioOut ( ofSoundBuffer& outBuffer )
{
	try {
		// Safety: Fill with silence first
		for ( size_t i = 0; i < outBuffer.getNumFrames ( ); i++ )
		{
			outBuffer.getSample ( i, 0 ) = 0.0f;
			outBuffer.getSample ( i, 1 ) = 0.0f;
		}
		
	
	// Safety check for required data structures
	if ( !mRawView || !mRawView->GetAudioData ( ) || !mRawView->GetDataset ( ) )
	{
		// Fill buffer with silence if data not available
		for ( size_t i = 0; i < outBuffer.getNumFrames ( ); i++ )
		{
			outBuffer.getSample ( i, 0 ) = 0.0f;
			outBuffer.getSample ( i, 1 ) = 0.0f;
		}
		return;
	}
	
	// Update morphing parameters if needed
	if ( mMorphParametersChanged )
	{
		UpdateMorphParameters ( );
	}
	
	if ( bResetFlag )
	{
		for ( size_t i = 0; i < mPlayheads.size ( ); i++ )
		{
			mPlayheads.erase ( mPlayheads.begin ( ) + i );
			i--;
		}

		bResetFlag = false;
	}

	std::vector<size_t> playheadsToKillThisBuffer;

	// -------------------------------------------------------------------------
	// NEW: Transfer queued playheads from the main thread to the audio thread
	//      and fetch any playheads that need to be killed. This restores the
	//      behaviour from the reference implementation that was inadvertently
	//      removed, and prevents the mNewPlayheads queue from filling up and
	//      rejecting new playhead requests.
	// -------------------------------------------------------------------------
	if ( mNewPlayheadMutex.try_lock ( ) )
	{
		std::lock_guard<std::mutex> lock ( mNewPlayheadMutex, std::adopt_lock );

		// Move newly created playheads into the active list
		while ( !mNewPlayheads.empty ( ) )
		{
			mPlayheads.push_back ( mNewPlayheads.front ( ) );
			mNewPlayheads.pop ( );
		}

		// Collect playheads that should be killed this buffer
		while ( !mPlayheadsToKill.empty ( ) )
		{
			playheadsToKillThisBuffer.push_back ( mPlayheadsToKill.front ( ) );
			mPlayheadsToKill.pop ( );
		}
	}

	// -------------------------------------------------------------------------
	// End of NEW section
	// -------------------------------------------------------------------------

	double crossoverJumpChance = (double)mCrossoverJumpChanceX1000 / 1000.0;

	for ( size_t playheadIndex = 0; playheadIndex < mPlayheads.size ( ); playheadIndex++ )
	{
		ofSoundBuffer playheadBuffer;
		playheadBuffer.setSampleRate ( mSampleRate );
		playheadBuffer.allocate ( outBuffer.getNumFrames ( ), 1 );
		
		for ( size_t i = 0; i < playheadBuffer.size ( ); i++ )
		{
			playheadBuffer.getSample ( i, 0 ) = 0.0;
		}

		size_t playheadBufferPosition = 0;
		// processing loop
		while ( true )
		{
			// crossfade jump
			// CrossfadeAudioSegment ( &playheadBuffer, &playheadBufferPosition, jumpOriginStartSample, jumpOriginEndSample, jumpOriginFile, &mPlayheads[playheadIndex], mCrossfadeSampleLength, true );
			if ( mPlayheads[playheadIndex].crossfading )
			{
				size_t crossfadeSamplesLeft = mPlayheads[playheadIndex].crossfadeSampleLength - mPlayheads[playheadIndex].crossfadeCurrentSample;
				size_t bufferSpace = playheadBuffer.getNumFrames ( ) - playheadBufferPosition;
				if ( crossfadeSamplesLeft > bufferSpace ) { crossfadeSamplesLeft = bufferSpace; }

				if ( mMorphMode && mAudioTransport && mAudioTransport->initialized ( ) )
				{
					// Use AudioTransport morphing
					ProcessMorphFrame ( &playheadBuffer, &playheadBufferPosition, &mPlayheads[playheadIndex], crossfadeSamplesLeft );
				}
				else
				{
					// Use traditional crossfade
					for ( size_t i = 0; i < crossfadeSamplesLeft; i++ )
					{
						float gain_A = cos ( (float)( mPlayheads[playheadIndex].crossfadeCurrentSample + i ) / (float)mPlayheads[playheadIndex].crossfadeSampleLength * 0.5 * M_PI );
						float gain_B = sin ( (float)( mPlayheads[playheadIndex].crossfadeCurrentSample + i ) / (float)mPlayheads[playheadIndex].crossfadeSampleLength * 0.5 * M_PI );

						float sample_A = 0.0f, sample_B = 0.0f;
						
						// Safely get sample A
						if ( mRawView && mRawView->GetAudioData ( ) &&
						     mPlayheads[playheadIndex].fileIndex < mRawView->GetAudioData ( )->raw.size() &&
						     mPlayheads[playheadIndex].fileIndex < mRawView->GetAudioData ( )->loaded.size() &&
						     mRawView->GetAudioData ( )->loaded[mPlayheads[playheadIndex].fileIndex] )
						{
							auto& bufferA = mRawView->GetAudioData ( )->raw[mPlayheads[playheadIndex].fileIndex];
							size_t idxA = mPlayheads[playheadIndex].sampleIndex + i;
							if ( idxA < bufferA.size() )
							{
								sample_A = bufferA.getSample ( idxA, 0 );
							}
						}
						
						// Safely get sample B
						if ( mRawView && mRawView->GetAudioData ( ) &&
						     mPlayheads[playheadIndex].jumpFileIndex < mRawView->GetAudioData ( )->raw.size() &&
						     mPlayheads[playheadIndex].jumpFileIndex < mRawView->GetAudioData ( )->loaded.size() &&
						     mRawView->GetAudioData ( )->loaded[mPlayheads[playheadIndex].jumpFileIndex] )
						{
							auto& bufferB = mRawView->GetAudioData ( )->raw[mPlayheads[playheadIndex].jumpFileIndex];
							size_t idxB = mPlayheads[playheadIndex].jumpSampleIndex + i;
							if ( idxB < bufferB.size() )
							{
								sample_B = bufferB.getSample ( idxB, 0 );
							}
						}
						
						playheadBuffer.getSample ( playheadBufferPosition + i, 0 ) = sample_A * gain_A + sample_B * gain_B;
					}

					mPlayheads[playheadIndex].crossfadeCurrentSample += crossfadeSamplesLeft;
					mPlayheads[playheadIndex].sampleIndex += crossfadeSamplesLeft;
					mPlayheads[playheadIndex].jumpSampleIndex += crossfadeSamplesLeft;
					playheadBufferPosition += crossfadeSamplesLeft;
				}

				if ( mPlayheads[playheadIndex].crossfadeCurrentSample >= mPlayheads[playheadIndex].crossfadeSampleLength )
				{
					mPlayheads[playheadIndex].crossfading = false;
					mPlayheads[playheadIndex].fileIndex = mPlayheads[playheadIndex].jumpFileIndex;
					mPlayheads[playheadIndex].sampleIndex = mPlayheads[playheadIndex].jumpSampleIndex;
					CalculateTriggerPoints ( mPlayheads[playheadIndex] );
				}
				else
				{
					break;
				}
			}

			// remove trigger points that have been hit
			while ( mPlayheads[playheadIndex].triggerSamplePoints.size ( ) > 0 && mPlayheads[playheadIndex].sampleIndex >= mPlayheads[playheadIndex].triggerSamplePoints.front ( ) )
			{
				mPlayheads[playheadIndex].triggerSamplePoints.pop ( );
			}

			// if EOF: loop/kill
			if ( mPlayheads[playheadIndex].triggerSamplePoints.size ( ) == 0 )
			{
				if ( mLoopPlayheads )
				{
					mPlayheads[playheadIndex].sampleIndex = 0;
					CalculateTriggerPoints ( mPlayheads[playheadIndex] );
				}
				else
				{
					playheadsToKillThisBuffer.push_back ( mPlayheads[playheadIndex].playheadID );
					break;
				}
			}

			// exit loop - no more space in outbuffer, no more triggers hit
			if ( !mPlayheads[playheadIndex].triggerSamplePoints.empty() && 
			     ( playheadBuffer.size ( ) - playheadBufferPosition ) < ( mPlayheads[playheadIndex].triggerSamplePoints.front ( ) - mPlayheads[playheadIndex].sampleIndex ) )
			{
				FillAudioSegment ( &playheadBuffer, &playheadBufferPosition, &mPlayheads[playheadIndex], true );
				break;
			}

			// fill audio up to the next trigger, already established that there is enough space in outBuffer
			FillAudioSegment ( &playheadBuffer, &playheadBufferPosition, &mPlayheads[playheadIndex], false );

			// after this point it is assumed that a new trigger has been reached, perform jump checks for this trigger
			
			int requiredSamples = mCrossfadeSampleLength;
			if ( mMorphMode )
			{
				// Convert duration in seconds to samples
				requiredSamples = mMorphTransitionDuration * mSampleRate;
			}
			// Safety check before accessing audio data
			if ( !mRawView || !mRawView->GetAudioData ( ) || 
			     mPlayheads[playheadIndex].fileIndex >= mRawView->GetAudioData ( )->raw.size() ||
			     mPlayheads[playheadIndex].sampleIndex + requiredSamples >= mRawView->GetAudioData ( )->raw[mPlayheads[playheadIndex].fileIndex].getNumFrames ( ) ) 
			{ 
				static int safetyFailCount = 0;
				if ( ++safetyFailCount % 100 == 0 )
				{
					ofLogWarning ( "AudioPlayback" ) << "Safety check failed " << safetyFailCount << " times. "
					                                  << "fileIndex=" << mPlayheads[playheadIndex].fileIndex 
					                                  << ", sampleIndex=" << mPlayheads[playheadIndex].sampleIndex
					                                  << ", requiredSamples=" << requiredSamples;
				}
				continue; 
			}
			
			double randomValue = (double)rand ( ) / RAND_MAX;
			if ( randomValue > crossoverJumpChance ) { continue; }
			
			// Different jump logic for time analysis vs stats mode
			if ( mRawView->IsTimeAnalysis() && mTimeCorpusMutex.try_lock ( ) )
			{
				std::lock_guard<std::mutex> lock ( mTimeCorpusMutex, std::adopt_lock );
				
				// Safety check for time corpus bounds
				if ( mPlayheads[playheadIndex].fileIndex >= mTimeCorpus.size() )
				{
					continue;
				}

				size_t hopSize = mRawView->GetDataset ( )->analysisSettings.windowFFTSize / mRawView->GetDataset ( )->analysisSettings.hopFraction;
				size_t timePointIndex = mPlayheads[playheadIndex].sampleIndex / hopSize;
				
				// Check vertex bounds
				if ( timePointIndex >= mTimeCorpus[mPlayheads[playheadIndex].fileIndex].getNumVertices() )
				{
					continue;
				}
				
				glm::vec3 playheadPosition = mTimeCorpus[mPlayheads[playheadIndex].fileIndex].getVertex ( timePointIndex );
				Utils::PointFT nearestPoint;
				Utils::PointFT currentPoint; currentPoint.file = mPlayheads[playheadIndex].fileIndex; currentPoint.time = timePointIndex;

				if ( !mPointPicker || 
				     !mPointPicker->FindNearestToPosition ( playheadPosition, nearestPoint, currentPoint, 
															mMaxJumpDistanceSpaceX1000, mMaxJumpTargets, mJumpSameFileAllowed, 
															mJumpSameFileMinTimeDiff, requiredSamples, *mRawView->GetAudioData ( ), hopSize ) )
				{ 
					ofLogWarning ( "AudioPlayback" ) << "FindNearestToPosition failed or no PointPicker";
					continue; 
				}
				
				// Load the target file if not already loaded
				if ( !LoadAudioFile ( nearestPoint.file ) )
				{
					ofLogWarning ( "AudioPlayback" ) << "Failed to load target file " << nearestPoint.file << " for TIME mode jump";
					continue;
				}

				mPlayheads[playheadIndex].crossfading = true;
				mPlayheads[playheadIndex].jumpFileIndex = nearestPoint.file;
				mPlayheads[playheadIndex].jumpSampleIndex = nearestPoint.time * hopSize;
				mPlayheads[playheadIndex].crossfadeCurrentSample = 0;
				mPlayheads[playheadIndex].crossfadeSampleLength = requiredSamples;
			}
			else if ( !mRawView->IsTimeAnalysis() )
			{
				// Stats mode: randomly jump to another file
				
				size_t numFiles = mRawView->GetDataset()->fileList.size();
				if ( numFiles > 1 )
				{
					// Build list of ALL candidate files (not just loaded ones)
					std::vector<size_t> candidates;
					for ( size_t i = 0; i < numFiles; i++ )
					{
						// Skip current file if same file jumps not allowed
						if ( !mJumpSameFileAllowed && i == mPlayheads[playheadIndex].fileIndex )
							continue;
							
						// Add all files as candidates - we'll load on demand
						candidates.push_back ( i );
					}
					
					// Pick a random candidate
					if ( !candidates.empty() )
					{
						size_t targetIndex = candidates[rand() % candidates.size()];
						
						// Load the target file if not already loaded
						if ( !LoadAudioFile ( targetIndex ) )
						{
							ofLogWarning ( "AudioPlayback" ) << "Failed to load target file " << targetIndex;
							continue;
						}
						
						// For stats mode, jump to beginning of the target file
						mPlayheads[playheadIndex].crossfading = true;
						mPlayheads[playheadIndex].jumpFileIndex = targetIndex;
						mPlayheads[playheadIndex].jumpSampleIndex = 0;
						mPlayheads[playheadIndex].crossfadeCurrentSample = 0;
						mPlayheads[playheadIndex].crossfadeSampleLength = requiredSamples;
					}
					else
					{
						ofLogWarning ( "AudioPlayback" ) << "No valid jump candidates found";
					}
				}
				else
				{
					ofLogWarning ( "AudioPlayback" ) << "Only " << numFiles << " files available, need at least 2 for jumps";
				}
			}
		}
		
		// if playhead is marked for death, apply a fade out and remove from playheads
		{
			std::vector<size_t>::iterator it = std::find ( playheadsToKillThisBuffer.begin ( ), playheadsToKillThisBuffer.end ( ), mPlayheads[playheadIndex].playheadID );
			size_t killIndex = std::distance ( playheadsToKillThisBuffer.begin ( ), it );
			if ( it != playheadsToKillThisBuffer.end ( ) )
			{
				for ( size_t i = 0; i < playheadBuffer.size ( ); i++ )
				{
					float gain = cos ( (float)i / (float)playheadBuffer.size ( ) * 0.5 * M_PI );
					playheadBuffer.getSample ( i, 0 ) *= gain;
				}

				playheadsToKillThisBuffer.erase ( playheadsToKillThisBuffer.begin ( ) + killIndex );

				mPlayheads.erase ( mPlayheads.begin ( ) + playheadIndex );
				playheadIndex--;
			}
		}

		// processing done, add to outBuffer
		float maxSample = 0.0f;
		for ( size_t sampleIndex = 0; sampleIndex < outBuffer.getNumFrames ( ); sampleIndex++ )
		{
			float sample = playheadBuffer.getSample ( sampleIndex, 0 );
			maxSample = std::max(maxSample, std::abs(sample));
			outBuffer.getSample ( sampleIndex, 0 ) += sample;
			outBuffer.getSample ( sampleIndex, 1 ) += sample;
		}
		
	}

	if ( mVisualPlayheadUpdateMutex.try_lock ( ) )
	{
		std::lock_guard<std::mutex> lock ( mVisualPlayheadUpdateMutex, std::adopt_lock );

		if ( mVisualPlayheads.size ( ) > 0 ) { mVisualPlayheads.clear ( ); }

		for ( size_t i = 0; i < mPlayheads.size ( ); i++ )
		{
			mVisualPlayheads.push_back ( Utils::VisualPlayhead ( mPlayheads[i].playheadID, mPlayheads[i].fileIndex, mPlayheads[i].sampleIndex ) );
		}
	}

	mActivePlayheads = mPlayheads.size ( );
	
	} catch ( const std::exception& e ) {
		ofLogError ( "AudioPlayback::audioOut" ) << "Exception caught: " << e.what();
	} catch ( ... ) {
		ofLogError ( "AudioPlayback::audioOut" ) << "Unknown exception caught";
	}
}

void Explorer::AudioPlayback::FillAudioSegment ( ofSoundBuffer* outBuffer, size_t* outBufferPosition, Utils::AudioPlayhead* playhead, bool outBufferFull )
{
	// If no trigger points, fill the remaining buffer or to end of audio file
	size_t segmentLength;
	if ( playhead->triggerSamplePoints.empty() )
	{
		// Fill to end of audio file or buffer, whichever is smaller
		if ( mRawView && mRawView->GetAudioData() && 
		     playhead->fileIndex < mRawView->GetAudioData()->raw.size() )
		{
			size_t samplesLeft = mRawView->GetAudioData()->raw[playhead->fileIndex].size() - playhead->sampleIndex;
			segmentLength = std::min(samplesLeft, outBuffer->size() - *outBufferPosition);
		}
		else
		{
			segmentLength = outBuffer->size() - *outBufferPosition;
		}
	}
	else
	{
		segmentLength = playhead->triggerSamplePoints.front ( ) - playhead->sampleIndex;
	}

	if ( outBufferFull && segmentLength > (outBuffer->size ( ) - *outBufferPosition) ) // cut off early if outBuffer is full
	{
		segmentLength = outBuffer->size ( ) - *outBufferPosition;
	}

	for ( size_t i = 0; i < segmentLength; i++ )
	{
		float sample = 0.0f;
		
		// Safely get sample with bounds checking
		if ( mRawView && mRawView->GetAudioData ( ) &&
		     playhead->fileIndex < mRawView->GetAudioData ( )->raw.size() &&
		     playhead->fileIndex < mRawView->GetAudioData ( )->loaded.size() &&
		     mRawView->GetAudioData ( )->loaded[playhead->fileIndex] )
		{
			auto& buffer = mRawView->GetAudioData ( )->raw[playhead->fileIndex];
			size_t idx = playhead->sampleIndex + i;
			if ( idx < buffer.size() )
			{
				sample = buffer.getSample ( idx, 0 );
			}
			else
			{
				// If we've reached the end of the file, no sample
				sample = 0.0f;
			}
		}
		else
		{
			// If file not loaded or invalid index, output silence
			static int warningCount = 0;
			if ( ++warningCount % 1000 == 0 )
			{
				ofLogWarning ( "AudioPlayback" ) << "FillAudioSegment: File " << playhead->fileIndex << " not loaded or invalid";
			}
		}
		
		outBuffer->getSample ( *outBufferPosition + i, 0 ) = sample;
	}

	playhead->sampleIndex += segmentLength;
	*outBufferPosition += segmentLength;
}

void Explorer::AudioPlayback::CrossfadeAudioSegment ( ofSoundBuffer* outBuffer, size_t* outBufferPosition, size_t startSample_A, size_t endSample_A, size_t fileIndex_A, Utils::AudioPlayhead* playhead_B, size_t lengthSetting, bool outBufferFull )
{
	size_t originLength = endSample_A - startSample_A;
	size_t jumpLength = playhead_B->triggerSamplePoints.front ( ) - playhead_B->sampleIndex;
	size_t crossfadeLength = (originLength < jumpLength) ? originLength : jumpLength;
	crossfadeLength = (crossfadeLength < lengthSetting) ? crossfadeLength : lengthSetting;

	if ( outBufferFull && crossfadeLength > ( outBuffer->size ( ) - *outBufferPosition ) ) // cut off early if outBuffer is full
	{
		crossfadeLength = outBuffer->size ( ) - *outBufferPosition;
	}


	for ( size_t i = 0; i < crossfadeLength; i++ )
	{
		float gain_A = cos ( (float)i / (float)crossfadeLength * 0.5 * M_PI );
		float gain_B = sin ( (float)i / (float)crossfadeLength * 0.5 * M_PI );

		float sample_A = 0.0f, sample_B = 0.0f;
		
		// Safely get sample A
		if ( mRawView && mRawView->GetAudioData ( ) &&
		     fileIndex_A < mRawView->GetAudioData ( )->raw.size() &&
		     fileIndex_A < mRawView->GetAudioData ( )->loaded.size() &&
		     mRawView->GetAudioData ( )->loaded[fileIndex_A] )
		{
			auto& bufferA = mRawView->GetAudioData ( )->raw[fileIndex_A];
			size_t idxA = startSample_A + i;
			if ( idxA < bufferA.size() )
			{
				sample_A = bufferA.getSample ( idxA, 0 );
			}
		}
		
		// Safely get sample B
		if ( mRawView && mRawView->GetAudioData ( ) &&
		     playhead_B->fileIndex < mRawView->GetAudioData ( )->raw.size() &&
		     playhead_B->fileIndex < mRawView->GetAudioData ( )->loaded.size() &&
		     mRawView->GetAudioData ( )->loaded[playhead_B->fileIndex] )
		{
			auto& bufferB = mRawView->GetAudioData ( )->raw[playhead_B->fileIndex];
			size_t idxB = playhead_B->sampleIndex + i;
			if ( idxB < bufferB.size() )
			{
				sample_B = bufferB.getSample ( idxB, 0 );
			}
		}
		
		outBuffer->getSample ( *outBufferPosition + i, 0 ) = sample_A * gain_A + sample_B * gain_B;
	}

	playhead_B->sampleIndex += crossfadeLength;
	*outBufferPosition += crossfadeLength;
}

bool Explorer::AudioPlayback::CreatePlayhead ( size_t fileIndex, size_t sampleIndex )
{
	// Load audio file first (before locking mutex to avoid holding lock during file I/O)
	if ( !LoadAudioFile ( fileIndex ) )
	{
		return false;
	}

	Utils::AudioPlayhead newPlayhead ( playheadCounter, fileIndex, sampleIndex );
	playheadCounter++;

	CalculateTriggerPoints ( newPlayhead );

	// Lock mutex once for both size check and push
	{
		std::lock_guard<std::mutex> lock ( mNewPlayheadMutex );
		if ( mNewPlayheads.size ( ) > 3 )
		{
			ofLogWarning ( "AudioPlayback" ) << "Too many playheads queued already, failed to create new playhead";
			return false;
		}
		
		mNewPlayheads.push ( newPlayhead );
	}

	return true;
}

bool Explorer::AudioPlayback::KillPlayhead ( size_t playheadID )
{
	{
		std::lock_guard<std::mutex> lock ( mNewPlayheadMutex );
		mPlayheadsToKill.push ( playheadID );
	}

	return true;
}

std::vector<Utils::VisualPlayhead> Explorer::AudioPlayback::GetPlayheadInfo ( )
{
	std::lock_guard<std::mutex> lock ( mVisualPlayheadUpdateMutex );

	return mVisualPlayheads;
}

void Explorer::AudioPlayback::SetFlagReset ( )
{
	bResetFlag = true;
}

void Explorer::AudioPlayback::WaitForResetConfirm ( )
{
	if ( !bStreamStarted ) { return; }
	while ( bResetFlag )
	{
		if ( mActivePlayheads == 0 ) { return; }
	}
}

void Explorer::AudioPlayback::SetTimeCorpus ( const std::vector<ofMesh>& timeCorpus )
{
	std::lock_guard<std::mutex> lock ( mTimeCorpusMutex );

	mTimeCorpus = timeCorpus;
}

void Explorer::AudioPlayback::CalculateTriggerPoints ( Utils::AudioPlayhead& playhead )
{
	while ( !playhead.triggerSamplePoints.empty ( ) )
	{
		playhead.triggerSamplePoints.pop ( );
	}

	// Check if mRawView is valid
	if ( !mRawView || !mRawView->GetDataset ( ) || !mRawView->GetAudioData ( ) )
	{
		// Add a single trigger point at the end if we can't calculate properly
		playhead.triggerSamplePoints.push ( INT_MAX );
		return;
	}

	// Additional safety check for raw data
	if ( playhead.fileIndex >= mRawView->GetAudioData ( )->raw.size() )
	{
		playhead.triggerSamplePoints.push ( INT_MAX );
		return;
	}

	int triggerPointDistance = mRawView->GetDataset ( )->analysisSettings.windowFFTSize / mRawView->GetDataset ( )->analysisSettings.hopFraction;
	int currentTriggerPoint = triggerPointDistance;

	while ( currentTriggerPoint < mRawView->GetAudioData ( )->raw[playhead.fileIndex].size ( ) )
	{
		if ( currentTriggerPoint > playhead.sampleIndex )
		{
			playhead.triggerSamplePoints.push ( currentTriggerPoint );
		}
		currentTriggerPoint += triggerPointDistance;
	}

	playhead.triggerSamplePoints.push ( mRawView->GetAudioData ( )->raw[playhead.fileIndex].size ( ) - 1 );
}

void Explorer::AudioPlayback::SetMorphSTFTSize ( int size )
{
	std::lock_guard<std::mutex> lock ( mMorphParamMutex );
	mStagedSTFTSize = size;
	mMorphParametersChanged = true;
}

void Explorer::AudioPlayback::UpdateMorphParameters ( )
{
	std::lock_guard<std::mutex> lock ( mMorphParamMutex );
	if ( mStagedSTFTSize != mMorphSTFTSize )
	{
		mMorphSTFTSize = mStagedSTFTSize;
		// Reinitialize AudioTransport with new size
		if ( mAudioTransport )
		{
			mAudioTransport->init ( mMorphSTFTSize, mMorphSTFTSize, mMorphSTFTSize / 2 );
		}
	}
	mMorphParametersChanged = false;
}

void Explorer::AudioPlayback::ProcessMorphFrame ( ofSoundBuffer* outBuffer, size_t* outBufferPosition, 
                                                  Utils::AudioPlayhead* playhead, size_t crossfadeLength )
{
	// Ensure we have valid data
	if ( !mRawView || !mRawView->GetAudioData ( ) )
	{
		// Fall back to silence
		for ( size_t i = 0; i < crossfadeLength; i++ )
		{
			outBuffer->getSample ( *outBufferPosition + i, 0 ) = 0.0f;
		}
		playhead->crossfadeCurrentSample += crossfadeLength;
		playhead->sampleIndex += crossfadeLength;
		playhead->jumpSampleIndex += crossfadeLength;
		*outBufferPosition += crossfadeLength;
		return;
	}
	
	// Get the current interpolation weight based on crossfade progress
	float weight = (float)playhead->crossfadeCurrentSample / (float)playhead->crossfadeSampleLength;
	
	// Get STFT frame size
	size_t frameSize = mMorphSTFTSize;
	
	// Prepare input frames using fluid::RealVector
	fluid::RealVector frame1 ( frameSize );
	fluid::RealVector frame2 ( frameSize );
	
	// Get audio data references
	auto& audioData = mRawView->GetAudioData ( )->raw;
	
	// Bounds checking for file indices
	if ( playhead->fileIndex >= audioData.size() || playhead->jumpFileIndex >= audioData.size() )
	{
		// Fall back to silence if indices are out of bounds
		for ( size_t i = 0; i < crossfadeLength; i++ )
		{
			outBuffer->getSample ( *outBufferPosition + i, 0 ) = 0.0f;
		}
		playhead->crossfadeCurrentSample += crossfadeLength;
		playhead->sampleIndex += crossfadeLength;
		playhead->jumpSampleIndex += crossfadeLength;
		*outBufferPosition += crossfadeLength;
		return;
	}
	
	size_t audioSize1 = audioData[playhead->fileIndex].size ( );
	size_t audioSize2 = audioData[playhead->jumpFileIndex].size ( );
	
	// Fill frames from audio data with bounds checking
	for ( size_t i = 0; i < frameSize; i++ )
	{
		size_t idx1 = playhead->sampleIndex + i;
		size_t idx2 = playhead->jumpSampleIndex + i;
		
		// fluid::RealVector uses double precision
		frame1[i] = ( idx1 < audioSize1 ) ? 
			static_cast<double>(audioData[playhead->fileIndex].getSample ( idx1, 0 )) : 0.0;
		frame2[i] = ( idx2 < audioSize2 ) ? 
			static_cast<double>(audioData[playhead->jumpFileIndex].getSample ( idx2, 0 )) : 0.0;
	}
	
	// Process morph
	mAudioTransport->processFrame ( frame1, frame2, weight, mMorphOutputBuffer, mAllocator );
	
	// Copy morphed output with window normalization
	auto morphed = mMorphOutputBuffer.row ( 0 );
	auto window = mMorphOutputBuffer.row ( 1 );
	
	for ( size_t i = 0; i < crossfadeLength && i < frameSize; i++ )
	{
		float sample = morphed[i] / ( window[i] > 0 ? window[i] : 1.0 );
		outBuffer->getSample ( *outBufferPosition + i, 0 ) = sample;
	}
	
	// Update positions
	playhead->crossfadeCurrentSample += crossfadeLength;
	playhead->sampleIndex += crossfadeLength;
	playhead->jumpSampleIndex += crossfadeLength;
	*outBufferPosition += crossfadeLength;
}

bool Explorer::AudioPlayback::LoadAudioFile ( size_t fileIndex )
{
	// Check if we have valid data structures
	if ( !mRawView || !mRawView->GetAudioData ( ) || !mRawView->GetDataset ( ) )
	{
		ofLogError ( "AudioPlayback" ) << "LoadAudioFile: Invalid data structures";
		return false;
	}
	
	auto audioData = mRawView->GetAudioData ( );
	auto dataset = mRawView->GetDataset ( );
	
	// Ensure vectors are properly sized
	if ( audioData->loaded.size() < dataset->fileList.size() )
	{
		audioData->loaded.resize ( dataset->fileList.size(), false );
		audioData->raw.resize ( dataset->fileList.size() );
	}
	
	// Check bounds
	if ( fileIndex >= dataset->fileList.size() )
	{
		ofLogError ( "AudioPlayback" ) << "File index " << fileIndex << " out of bounds (max: " << dataset->fileList.size() - 1 << ")";
		return false;
	}
	
	// Check if already loaded
	if ( audioData->loaded[fileIndex] )
	{
		return true;
	}
	
	ofLogNotice ( "AudioPlayback" ) << "Loading audio file: " << dataset->fileList[fileIndex];
	
	// Use the AudioFileLoader class to load the file properly
	fluid::RealVector audioVector;
	if ( !Utils::AudioFileLoader::ReadAudioFile ( dataset->fileList[fileIndex], audioVector, mSampleRate ) )
	{
		ofLogError ( "AudioPlayback" ) << "Failed to load audio file: " << dataset->fileList[fileIndex];
		return false;
	}
	
	// Convert fluid::RealVector to ofSoundBuffer
	// AudioFileLoader returns mono audio at the target sample rate
	audioData->raw[fileIndex].allocate ( audioVector.size(), 1 );
	audioData->raw[fileIndex].setSampleRate( mSampleRate );
	
	// Copy the audio data
	std::vector<float> audioBuffer(audioVector.size());
	for ( size_t i = 0; i < audioVector.size(); i++ )
	{
		audioBuffer[i] = static_cast<float>(audioVector[i]);
	}
	audioData->raw[fileIndex].copyFrom(audioBuffer.data(), audioBuffer.size(), 1, mSampleRate);
	
	audioData->loaded[fileIndex] = true;
	ofLogNotice ( "AudioPlayback" ) << "Successfully loaded audio file";
	
	return true;
}

void Explorer::AudioPlayback::SetMorphTargets ( const std::vector<std::pair<size_t, size_t>>& targets,
                                               const AudioTransportN::BarycentricWeights& weights )
{
	if ( targets.size() != weights.size() ) {
		ofLogError ( "AudioPlayback" ) << "SetMorphTargets: targets and weights size mismatch";
		return;
	}
	
	if ( targets.size() > AudioTransportN::getMaxSources() ) {
		ofLogError ( "AudioPlayback" ) << "SetMorphTargets: too many targets, max is " 
			<< AudioTransportN::getMaxSources();
		return;
	}
	
	std::lock_guard<std::mutex> lock ( mMorphTargetsMutex );
	mMorphTargets = targets;
	mMorphWeights = weights;
	mMorphWeights.normalize(); // Ensure weights sum to 1.0
}

void Explorer::AudioPlayback::ClearMorphTargets ( )
{
	std::lock_guard<std::mutex> lock ( mMorphTargetsMutex );
	mMorphTargets.clear();
	mMorphWeights = AudioTransportN::BarycentricWeights();
}

void Explorer::AudioPlayback::ProcessNWayMorphFrame ( ofSoundBuffer* outBuffer, size_t* outBufferPosition,
                                                     Utils::AudioPlayhead* playhead, size_t crossfadeLength )
{
	// Ensure we have valid data
	if ( !mRawView || !mRawView->GetAudioData ( ) || !mAudioTransportN || !mAudioTransportN->initialized() ) {
		// Fall back to silence
		for ( size_t i = 0; i < crossfadeLength; i++ ) {
			outBuffer->getSample ( *outBufferPosition + i, 0 ) = 0.0f;
		}
		*outBufferPosition += crossfadeLength;
		return;
	}
	
	// Get current morph targets
	std::vector<std::pair<size_t, size_t>> targets;
	AudioTransportN::BarycentricWeights weights;
	{
		std::lock_guard<std::mutex> lock ( mMorphTargetsMutex );
		targets = mMorphTargets;
		weights = mMorphWeights;
	}
	
	// Validate we have enough targets
	if ( targets.size() < 2 ) {
		// Fall back to standard morphing if we have exactly 2 targets
		if ( targets.size() == 2 && mAudioTransport && mAudioTransport->initialized() ) {
			// Use the second target's weight as the interpolation parameter
			float weight2 = weights[1];
			
			// TODO: Implement conversion to standard morph frame processing
			// For now, fall back to silence
			for ( size_t i = 0; i < crossfadeLength; i++ ) {
				outBuffer->getSample ( *outBufferPosition + i, 0 ) = 0.0f;
			}
			*outBufferPosition += crossfadeLength;
			return;
		}
		
		// Otherwise output silence
		for ( size_t i = 0; i < crossfadeLength; i++ ) {
			outBuffer->getSample ( *outBufferPosition + i, 0 ) = 0.0f;
		}
		*outBufferPosition += crossfadeLength;
		return;
	}
	
	// Get STFT frame size
	size_t frameSize = mMorphSTFTSize;
	
	// Prepare input frames
	std::vector<fluid::RealVector> frames;
	frames.reserve(targets.size());
	
	// Get audio data references
	auto& audioData = mRawView->GetAudioData ( )->raw;
	
	// Collect frames from each target
	for ( const auto& target : targets ) {
		size_t fileIndex = target.first;
		size_t sampleIndex = target.second;
		
		// Bounds checking
		if ( fileIndex >= audioData.size() ) {
			// Add silence frame
			frames.emplace_back(frameSize);
			for ( size_t i = 0; i < frameSize; i++ ) {
				frames.back()[i] = 0.0;
			}
			continue;
		}
		
		size_t audioSize = audioData[fileIndex].size();
		frames.emplace_back(frameSize);
		
		// Fill frame from audio data
		for ( size_t i = 0; i < frameSize; i++ ) {
			size_t idx = sampleIndex + i;
			frames.back()[i] = ( idx < audioSize ) ? 
				static_cast<double>(audioData[fileIndex].getSample ( idx, 0 )) : 0.0;
		}
	}
	
	// Convert frames to views
	std::vector<fluid::RealVectorView> frameViews;
	frameViews.reserve(frames.size());
	for ( auto& frame : frames ) {
		frameViews.push_back(frame);
	}
	
	// Process N-way morph
	mAudioTransportN->processFrameN(frameViews, weights, mMorphOutputBuffer);
	
	// Copy morphed output with window normalization
	auto morphed = mMorphOutputBuffer.row ( 0 );
	auto window = mMorphOutputBuffer.row ( 1 );
	
	for ( size_t i = 0; i < crossfadeLength && i < frameSize; i++ ) {
		float sample = morphed[i] / ( window[i] > 0 ? window[i] : 1.0 );
		outBuffer->getSample ( *outBufferPosition + i, 0 ) = sample;
	}
	
	// Update position
	*outBufferPosition += crossfadeLength;
	
	// Update playhead sample indices for all targets
	// This is a simplified approach - in practice you might want more sophisticated
	// playhead management for N-way morphing
	if ( playhead ) {
		playhead->sampleIndex += crossfadeLength;
	}
}