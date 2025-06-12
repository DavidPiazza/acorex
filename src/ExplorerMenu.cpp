#include "./ExplorerMenu.h"
#include <ofUtils.h>
#include <cstring>

using namespace Acorex;

void ExplorerMenu::Initialise ( bool HiDpi )
{
	// DPI -----------------------------------------
    {
		if ( HiDpi ) { mLayout.enableHiDpi ( ); }
		else { mLayout.disableHiDpi ( ); }
	}
	
	// Pointer Sharing -----------------------------
	{
		if ( !bViewPointerShared )
		{
			mRawView = std::make_shared<Explorer::RawView> ( );
			mLiveView.SetRawView ( mRawView );
			bViewPointerShared = true;
		}
	}

	// Clear --------------------------------------
	{
		RemoveListeners ( );

		mMainPanel.clear ( );
		mMorphPanel.clear ( );
		mCrossoverPanel.clear ( );
		mAudioDevicePanel.clear ( );

		mLiveView.Initialise ( );
	}

	// Variables ----------------------------------
	{
		mLastUpdateTime = 0;
		mSlowUpdateInterval = 100;

		mOpenCorpusButtonClickTime = 0;
		mOpenCorpusButtonTimeout = 3000;
	}

	// States ------------------------------------
	{
		bIsCorpusOpen = false;
		bOpenCorpusDrawWarning = false;
	}

	// Main Panel --------------------------------
	{
		int dropdownScrollSpeed = 32;

		mMainPanel.setup ( );

		mMainPanel.add ( mCorpusNameLabel.setup ( "", bInitialiseShouldLoad ? mRawView->GetCorpusName ( ) : "No Corpus Loaded" ) );
		mMainPanel.add ( mOpenCorpusButton.setup ( "Open Corpus" ) );

		mCorpusNameLabel.setBackgroundColor ( mColors.interfaceBackgroundColor );
		mOpenCorpusButton.setBackgroundColor ( mColors.interfaceBackgroundColor );

		mDimensionDropdownX.reset ( );
		mDimensionDropdownY.reset ( );
		mDimensionDropdownZ.reset ( );
		mDimensionDropdownColor.reset ( );

		mDimensionDropdownX = make_unique<ofxDropdown> ( "X Dimension", dropdownScrollSpeed );
		mDimensionDropdownY = make_unique<ofxDropdown> ( "Y Dimension", dropdownScrollSpeed );
		mDimensionDropdownZ = make_unique<ofxDropdown> ( "Z Dimension", dropdownScrollSpeed );
		mDimensionDropdownColor = make_unique<ofxDropdown> ( "Color Dimension", dropdownScrollSpeed );

		if ( bInitialiseShouldLoad )
		{
			mDimensionDropdownX->add ( "None" );
			mDimensionDropdownY->add ( "None" );
			mDimensionDropdownZ->add ( "None" );
			mDimensionDropdownColor->add ( "None" );

			for ( auto& dimension : mRawView->GetDimensions ( ) )
			{
				mDimensionDropdownX->add ( dimension );
				mDimensionDropdownY->add ( dimension );
				mDimensionDropdownZ->add ( dimension );
				mDimensionDropdownColor->add ( dimension );
			}

			bool needStatisticDropdowns = !mRawView->IsTimeAnalysis ( ) && !mRawView->IsReduction ( );
			
			// Note: Transport analysis dimensions are not yet supported in the UI
			if ( mRawView->IsTransportAnalysis ( ) && mRawView->GetDimensions ( ).size ( ) == 0 )
			{
				ofLogWarning ( "ExplorerMenu" ) << "Transport analysis detected but no Transport dimensions available in UI yet";
			}

			mMainPanel.add ( mDimensionDropdownX.get ( ) );
			mMainPanel.add ( mDimensionDropdownY.get ( ) );
			mMainPanel.add ( mDimensionDropdownZ.get ( ) );
			mMainPanel.add ( mDimensionDropdownColor.get ( ) );
		}

		mDimensionDropdownX->disableMultipleSelection ( );
		mDimensionDropdownY->disableMultipleSelection ( );
		mDimensionDropdownZ->disableMultipleSelection ( );
		mDimensionDropdownColor->disableMultipleSelection ( );

		mDimensionDropdownX->enableCollapseOnSelection ( );
		mDimensionDropdownY->enableCollapseOnSelection ( );
		mDimensionDropdownZ->enableCollapseOnSelection ( );
		mDimensionDropdownColor->enableCollapseOnSelection ( );

		mDimensionDropdownX->setDropDownPosition ( ofxDropdown::DD_LEFT );
		mDimensionDropdownY->setDropDownPosition ( ofxDropdown::DD_LEFT );
		mDimensionDropdownZ->setDropDownPosition ( ofxDropdown::DD_LEFT );
		mDimensionDropdownColor->setDropDownPosition ( ofxDropdown::DD_LEFT );

		mDimensionDropdownX->setBackgroundColor ( mColors.interfaceBackgroundColor );
		mDimensionDropdownY->setBackgroundColor ( mColors.interfaceBackgroundColor );
		mDimensionDropdownZ->setBackgroundColor ( mColors.interfaceBackgroundColor );
		mDimensionDropdownColor->setBackgroundColor ( mColors.interfaceBackgroundColor );

		if ( bInitialiseShouldLoad )
		{
			mMainPanel.add ( mColorSpectrumSwitcher.setup ( "Color Spectrum: Red<->Blue", false ) );
			mColorSpectrumSwitcher.setBackgroundColor ( mColors.interfaceBackgroundColor );
		}

		mMainPanel.setPosition ( ofGetWidth ( ) - mLayout.explorePanelWidth, mLayout.explorePanelOriginY );
		mMainPanel.setWidthElements ( mLayout.explorePanelWidth );
		mMainPanel.disableHeader ( );
	}
	
	// Morph Panel Setup ------------------------
	{
		mMorphPanel.setup("Morph Controls");
		mMorphPanel.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		// Basic morph controls
		mMorphModeToggle.setup("Morph Mode", false); // Disabled due to Eigen assertion
		mMorphModeToggle.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mContinuousMorphModeToggle.setup("Continuous Morph Mode", false);
		mContinuousMorphModeToggle.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mKDTreeNeighborsField.setup("KD-Tree Neighbors", 3, 1, 20);
		mKDTreeNeighborsField.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mTransitionDurationSlider.setup("Transition Duration", 0.5, 0.1, 5.0);
		mTransitionDurationSlider.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mSTFTSizeSlider.setup("STFT Size", 1024, 512, 4096);
		mSTFTSizeSlider.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		// Interpolation and loop controls
		mInterpolationPositionSlider.setup("Morph Position", 0.0, 0.0, 1.0);
		mInterpolationPositionSlider.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mLoopModeToggle.setup("Loop Morph", false);
		mLoopModeToggle.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mLoadingLabel.setup("Status", "Ready");
		mLoadingLabel.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		// Add controls to panel
		mMorphPanel.add(&mMorphModeToggle);
		mMorphPanel.add(&mContinuousMorphModeToggle);
		mMorphPanel.add(&mKDTreeNeighborsField);
		mMorphPanel.add(&mTransitionDurationSlider);
		mMorphPanel.add(&mSTFTSizeSlider);
		mMorphPanel.add(&mInterpolationPositionSlider);
		mMorphPanel.add(&mLoopModeToggle);
		mMorphPanel.add(&mLoadingLabel);
		
		// Position the morph panel below the main panel
		mMorphPanel.setPosition(ofGetWidth() - mLayout.explorePanelWidth, 
			mMainPanel.getPosition().y + mMainPanel.getHeight() + 20);
		mMorphPanel.setWidthElements(mLayout.explorePanelWidth);
	}
	
	// Crossover Jump Panel Setup ----------------
	{
		mCrossoverPanel.setup("Crossover Jump Controls");
		mCrossoverPanel.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		// Initialize controls with default values
		mLoopPlayheadsToggle.setup("Loop Playheads", false);
		mLoopPlayheadsToggle.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mJumpSameFileAllowedToggle.setup("Allow Same File Jumps", false);
		mJumpSameFileAllowedToggle.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mJumpSameFileMinTimeDiffSlider.setup("Min Time Diff (s)", 2, 0, 10);
		mJumpSameFileMinTimeDiffSlider.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mCrossoverJumpChanceSlider.setup("Jump Chance (%)", 5.0, 0.0, 100.0);
		mCrossoverJumpChanceSlider.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mCrossfadeMaxSampleLengthSlider.setup("Crossfade Length", 256, 64, 2048);
		mCrossfadeMaxSampleLengthSlider.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mMaxJumpDistanceSpaceSlider.setup("Max Jump Distance (%)", 5.0, 0.0, 100.0);
		mMaxJumpDistanceSpaceSlider.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		mMaxJumpTargetsSlider.setup("Max Jump Targets", 5, 1, 20);
		mMaxJumpTargetsSlider.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		// Add controls to panel
		mCrossoverPanel.add(&mLoopPlayheadsToggle);
		mCrossoverPanel.add(&mJumpSameFileAllowedToggle);
		mCrossoverPanel.add(&mJumpSameFileMinTimeDiffSlider);
		mCrossoverPanel.add(&mCrossoverJumpChanceSlider);
		mCrossoverPanel.add(&mCrossfadeMaxSampleLengthSlider);
		mCrossoverPanel.add(&mMaxJumpDistanceSpaceSlider);
		mCrossoverPanel.add(&mMaxJumpTargetsSlider);
		
		// Position below morph panel
		mCrossoverPanel.setPosition(ofGetWidth() - mLayout.explorePanelWidth,
			mMorphPanel.getPosition().y + mMorphPanel.getHeight() + 20);
		mCrossoverPanel.setWidthElements(mLayout.explorePanelWidth);
		
		// Initialize AudioPlayback with default values
		if ( mLiveView.GetAudioPlayback ( ) )
		{
			mLiveView.GetAudioPlayback ( )->SetLoopPlayheads ( false );
			mLiveView.GetAudioPlayback ( )->SetJumpSameFileAllowed ( false );
			mLiveView.GetAudioPlayback ( )->SetJumpSameFileMinTimeDiff ( 2 );
			mLiveView.GetAudioPlayback ( )->SetCrossoverJumpChance ( 50 ); // 5% = 50 per thousand
			mLiveView.GetAudioPlayback ( )->SetCrossfadeSampleLength ( 256 );
			mLiveView.GetAudioPlayback ( )->SetMaxJumpDistanceSpace ( 50 ); // 5% = 50 per thousand
			mLiveView.GetAudioPlayback ( )->SetMaxJumpTargets ( 5 );
		}
	}
	
	// Audio Device Panel Setup -------------------
	{
		int dropdownScrollSpeed = 32;
		
		mAudioDevicePanel.setup("Audio Device Settings");
		mAudioDevicePanel.setBackgroundColor(mColors.interfaceBackgroundColor);
		
		// Get audio devices
		ofSoundStream tempStream;
		mAudioDevices = tempStream.getDeviceList ( ofSoundDevice::Api::DEFAULT );
		
		// Buffer size dropdown
		mBufferSizeDropdown.reset();
		mBufferSizeDropdown = make_unique<ofxDropdown>("Buffer Size", dropdownScrollSpeed);
		
		for (int i = 0; i < mBufferSizes.size(); i++) {
			mBufferSizeDropdown->add(std::to_string(mBufferSizes[i]));
		}
		mBufferSizeDropdown->setSelectedValueByName("256", true);
		mBufferSizeDropdown->disableMultipleSelection();
		mBufferSizeDropdown->enableCollapseOnSelection();
		mBufferSizeDropdown->setDropDownPosition(ofxDropdown::DD_LEFT);
		mBufferSizeDropdown->setBackgroundColor(mColors.interfaceBackgroundColor);
		
		// Output device dropdown
		mOutDeviceDropdown.reset();
		mOutDeviceDropdown = make_unique<ofxDropdown>("Output Device", dropdownScrollSpeed);
		
		for (int i = 0; i < mAudioDevices.size(); i++) {
			if (mAudioDevices[i].outputChannels > 0) {
				mOutDeviceDropdown->add(mAudioDevices[i].name);
			}
		}
		
		if (mAudioDevices.size() > 0) {
			// Select default device
			for (int i = 0; i < mAudioDevices.size(); i++) {
				if (mAudioDevices[i].isDefaultOutput) {
					mOutDeviceDropdown->setSelectedValueByName(mAudioDevices[i].name, true);
					mCurrentDeviceIndex = i;
					break;
				}
			}
		}
		
		mOutDeviceDropdown->disableMultipleSelection();
		mOutDeviceDropdown->enableCollapseOnSelection();
		mOutDeviceDropdown->setDropDownPosition(ofxDropdown::DD_LEFT);
		mOutDeviceDropdown->setBackgroundColor(mColors.interfaceBackgroundColor);
		
		// Add controls to panel
		mAudioDevicePanel.add(mBufferSizeDropdown.get());
		mAudioDevicePanel.add(mOutDeviceDropdown.get());
		
		// Position below crossover panel
		mAudioDevicePanel.setPosition(ofGetWidth() - mLayout.explorePanelWidth,
			mCrossoverPanel.getPosition().y + mCrossoverPanel.getHeight() + 20);
		mAudioDevicePanel.setWidthElements(mLayout.explorePanelWidth);
	}

	// Listeners --------------------------------
	{
		mOpenCorpusButton.addListener ( this, &ExplorerMenu::OpenCorpus );
		mDimensionDropdownX->addListener ( this, &ExplorerMenu::SwapDimensionX );
		mDimensionDropdownY->addListener ( this, &ExplorerMenu::SwapDimensionY );
		mDimensionDropdownZ->addListener ( this, &ExplorerMenu::SwapDimensionZ );
		mDimensionDropdownColor->addListener ( this, &ExplorerMenu::SwapDimensionColor );
		mColorSpectrumSwitcher.addListener ( this, &ExplorerMenu::SwitchColorSpectrum );
		
		// Morph Mode Listeners
		mMorphModeToggle.addListener ( this, &ExplorerMenu::onMorphModeChanged );
		mContinuousMorphModeToggle.addListener ( this, &ExplorerMenu::onContinuousMorphModeChanged );
		mKDTreeNeighborsField.addListener ( this, &ExplorerMenu::onKDTreeNeighborsChanged );
		mTransitionDurationSlider.addListener ( this, &ExplorerMenu::onTransitionDurationChanged );
		mSTFTSizeSlider.addListener ( this, &ExplorerMenu::onSTFTSizeChanged );
		mInterpolationPositionSlider.addListener ( this, &ExplorerMenu::onInterpolationPositionChanged );
		mLoopModeToggle.addListener ( this, &ExplorerMenu::onLoopModeChanged );
		
		// Crossover Jump Listeners
		mLoopPlayheadsToggle.addListener ( this, &ExplorerMenu::onLoopPlayheadsChanged );
		mJumpSameFileAllowedToggle.addListener ( this, &ExplorerMenu::onJumpSameFileAllowedChanged );
		mJumpSameFileMinTimeDiffSlider.addListener ( this, &ExplorerMenu::onJumpSameFileMinTimeDiffChanged );
		mCrossoverJumpChanceSlider.addListener ( this, &ExplorerMenu::onCrossoverJumpChanceChanged );
		mCrossfadeMaxSampleLengthSlider.addListener ( this, &ExplorerMenu::onCrossfadeMaxSampleLengthChanged );
		mMaxJumpDistanceSpaceSlider.addListener ( this, &ExplorerMenu::onMaxJumpDistanceSpaceChanged );
		mMaxJumpTargetsSlider.addListener ( this, &ExplorerMenu::onMaxJumpTargetsChanged );
		
		// Audio Device Listeners
		mBufferSizeDropdown->addListener ( this, &ExplorerMenu::onBufferSizeChanged );
		mOutDeviceDropdown->addListener ( this, &ExplorerMenu::onOutDeviceChanged );
		
		bListenersAdded = true;
	}

	bInitialiseShouldLoad = false;
}

void ExplorerMenu::Show ( )
{
	bDraw = true;
}

// hides menu without resetting
void ExplorerMenu::Hide ( )
{
	bDraw = false;
}

void ExplorerMenu::Draw ( )
{
	if ( !bDraw ) { return; }

	mLiveView.Draw ( );

	// call pointpicker draw

	mMainPanel.draw ( );
	mMorphPanel.draw ( );
	mCrossoverPanel.draw ( );
	mAudioDevicePanel.draw ( );
}

void ExplorerMenu::Update ( )
{
	mLiveView.Update ( );
	
	// TODO: Update morph UI state when morph mode is implemented in AudioPlayback
	/*if (mMorphModeToggle && bMorphMode) {
		// Update interpolation slider to reflect actual morph state
		if (mLiveView.isInSingleSoundMode()) {
			// In single sound mode, show interpolation at 0
			mInterpolationPositionSlider = 0.0f;
		} else {
			// In morph mode, show actual interpolation weight
			mInterpolationPositionSlider = mLiveView.getCurrentInterpolationWeight();
		}
		
		// Update loading status
		if (mLiveView.isLoadingTarget()) {
			mLoadingLabel = "Loading...";
		} else {
			mLoadingLabel = "Ready";
		}
	}*/

	if ( ofGetElapsedTimeMillis ( ) - mLastUpdateTime > mSlowUpdateInterval )
	{
		mLastUpdateTime = ofGetElapsedTimeMillis ( );
		SlowUpdate ( );
	}
}

void ExplorerMenu::SlowUpdate ( )
{
	if ( bOpenCorpusDrawWarning && ofGetElapsedTimeMillis ( ) - mOpenCorpusButtonClickTime > mOpenCorpusButtonTimeout )
	{
		bOpenCorpusDrawWarning = false;
		mOpenCorpusButton.setName ( "Open Corpus" );
	}

	mLiveView.SlowUpdate ( );
}

void ExplorerMenu::Exit ( )
{
	RemoveListeners ( );
	mLiveView.Exit ( );
}

void ExplorerMenu::RemoveListeners ( )
{
	if ( !bListenersAdded ) { return; }
	mOpenCorpusButton.removeListener ( this, &ExplorerMenu::OpenCorpus );
	mDimensionDropdownX->removeListener ( this, &ExplorerMenu::SwapDimensionX );
	mDimensionDropdownY->removeListener ( this, &ExplorerMenu::SwapDimensionY );
	mDimensionDropdownZ->removeListener ( this, &ExplorerMenu::SwapDimensionZ );
	mDimensionDropdownColor->removeListener ( this, &ExplorerMenu::SwapDimensionColor );
	mColorSpectrumSwitcher.removeListener ( this, &ExplorerMenu::SwitchColorSpectrum );
	
	// Remove Morph Mode Listeners
	mMorphModeToggle.removeListener ( this, &ExplorerMenu::onMorphModeChanged );
	mTransitionDurationSlider.removeListener ( this, &ExplorerMenu::onTransitionDurationChanged );
	mSTFTSizeSlider.removeListener ( this, &ExplorerMenu::onSTFTSizeChanged );
	mInterpolationPositionSlider.removeListener ( this, &ExplorerMenu::onInterpolationPositionChanged );
	mLoopModeToggle.removeListener ( this, &ExplorerMenu::onLoopModeChanged );
	
	// Remove Crossover Jump Listeners
	mLoopPlayheadsToggle.removeListener ( this, &ExplorerMenu::onLoopPlayheadsChanged );
	mJumpSameFileAllowedToggle.removeListener ( this, &ExplorerMenu::onJumpSameFileAllowedChanged );
	mJumpSameFileMinTimeDiffSlider.removeListener ( this, &ExplorerMenu::onJumpSameFileMinTimeDiffChanged );
	mCrossoverJumpChanceSlider.removeListener ( this, &ExplorerMenu::onCrossoverJumpChanceChanged );
	mCrossfadeMaxSampleLengthSlider.removeListener ( this, &ExplorerMenu::onCrossfadeMaxSampleLengthChanged );
	mMaxJumpDistanceSpaceSlider.removeListener ( this, &ExplorerMenu::onMaxJumpDistanceSpaceChanged );
	mMaxJumpTargetsSlider.removeListener ( this, &ExplorerMenu::onMaxJumpTargetsChanged );
	
	// Remove Audio Device Listeners
	mBufferSizeDropdown->removeListener ( this, &ExplorerMenu::onBufferSizeChanged );
	mOutDeviceDropdown->removeListener ( this, &ExplorerMenu::onOutDeviceChanged );
	
	bListenersAdded = false;
}

// Main Functions ------------------------------

void ExplorerMenu::OpenCorpus ( )
{
	bBlockDimensionFilling = true;

	if ( bIsCorpusOpen && !bOpenCorpusDrawWarning )
	{
		bOpenCorpusDrawWarning = true;
		mOpenCorpusButtonClickTime = ofGetElapsedTimeMillis ( );
		mOpenCorpusButton.setName ( "!! Close Current? !!" );
		return;
	}

	bOpenCorpusDrawWarning = false;
	mOpenCorpusButton.setName ( "Open Corpus" );

	bool success = mRawView->LoadCorpus ( );
	if ( !success ) { return; }
	
	// Check what types of analysis the corpus contains
	bool hasTimeAnalysis = mRawView->IsTimeAnalysis ( );
	bool hasTransportAnalysis = mRawView->IsTransportAnalysis ( );
	bool hasTransportData = mRawView->HasTransportData ( );
	
	if ( hasTimeAnalysis )
	{
		ofLogNotice ( "ExplorerMenu" ) << "Corpus has Time analysis";
	}
	if ( hasTransportAnalysis )
	{
		ofLogNotice ( "ExplorerMenu" ) << "Corpus has Transport analysis enabled";
		if ( hasTransportData )
		{
			auto transportData = mRawView->GetTransportData ( );
			ofLogNotice ( "ExplorerMenu" ) << "Transport data present for " << transportData->fileCount() << " files";
		}
		else
		{
			ofLogWarning ( "ExplorerMenu" ) << "Transport analysis enabled but no Transport data found";
		}
	}
	
	// Enable/disable Continuous Morph Mode based on analysis types
	bool canUseContinuousMorph = hasTimeAnalysis && hasTransportAnalysis && hasTransportData;
	mContinuousMorphModeToggle.setBackgroundColor(canUseContinuousMorph ? 
		mColors.interfaceBackgroundColor : ofColor(60, 60, 60));
	mContinuousMorphModeToggle.setTextColor(canUseContinuousMorph ? 
		ofColor(255) : ofColor(120));
	
	// Also configure the KD-tree neighbors field
	mKDTreeNeighborsField.setBackgroundColor(canUseContinuousMorph ? 
		mColors.interfaceBackgroundColor : ofColor(60, 60, 60));
	mKDTreeNeighborsField.setTextColor(canUseContinuousMorph ? 
		ofColor(255) : ofColor(120));
	
	if (!canUseContinuousMorph)
	{
		mContinuousMorphModeToggle = false; // Ensure it's off if not available
		ofLogNotice("ExplorerMenu") << "Continuous Morph Mode disabled - requires both Time and Transport analysis";
	}
	else
	{
		ofLogNotice("ExplorerMenu") << "Continuous Morph Mode available";
	}
	
	bInitialiseShouldLoad = true;
	Initialise ( mLayout.HiDpi );

	mLiveView.CreatePoints ( );

	// set default dropdown values
	std::string xDimension = "None", yDimension = "None", zDimension = "None", colorDimension = "None";
	{
		int dimensionCount = mRawView->GetDimensions ( ).size ( );
		if ( mRawView->IsTimeAnalysis ( ) || mRawView->IsReduction ( ) )
		{
			xDimension = mRawView->GetDimensions ( ).size ( ) > 0 ? mRawView->GetDimensions ( )[0] : "None";
			yDimension = mRawView->GetDimensions ( ).size ( ) > 1 ? mRawView->GetDimensions ( )[1] : "None";
			zDimension = mRawView->GetDimensions ( ).size ( ) > 2 ? mRawView->GetDimensions ( )[2] : "None";
			colorDimension = mRawView->GetDimensions ( ).size ( ) > 3 ? mRawView->GetDimensions ( )[3] : "None";
		}
		else
		{
			xDimension = mRawView->GetDimensions ( ).size ( ) > 0 ? mRawView->GetDimensions ( )[0] : "None";
			yDimension = mRawView->GetDimensions ( ).size ( ) > 7 ? mRawView->GetDimensions ( )[7] : "None";
			zDimension = mRawView->GetDimensions ( ).size ( ) > 14 ? mRawView->GetDimensions ( )[14] : "None";
			colorDimension = mRawView->GetDimensions ( ).size ( ) > 21 ? mRawView->GetDimensions ( )[21] : "None";
		}
		
	}

	mDimensionDropdownX->setSelectedValueByName ( xDimension, true );
	mDimensionDropdownY->setSelectedValueByName ( yDimension, true );
	mDimensionDropdownZ->setSelectedValueByName ( zDimension, true );
	mDimensionDropdownColor->setSelectedValueByName ( colorDimension, true );

	bBlockDimensionFilling = false;

	SwapDimension ( xDimension, Utils::Axis::X );
	SwapDimension ( yDimension, Utils::Axis::Y );
	SwapDimension ( zDimension, Utils::Axis::Z );

	bIsCorpusOpen = true;

	SwapDimension ( colorDimension, Utils::Axis::COLOR );
}

void ExplorerMenu::SwapDimension ( string dimension, Utils::Axis axis )
{
	if ( bBlockDimensionFilling ) { return; }

	if ( dimension == "None" )					{ mLiveView.FillDimensionNone ( axis ); }
	else
	{
		int dimensionIndex = GetDimensionIndex ( dimension );
		if ( dimensionIndex == -1 ) { return; }

		if ( mRawView->IsTimeAnalysis ( ) )		{ mLiveView.FillDimensionTime ( dimensionIndex, axis ); }
		else if ( !mRawView->IsReduction ( ) )	{ mLiveView.FillDimensionStats ( dimensionIndex, axis ); }
		else									{ mLiveView.FillDimensionStatsReduced ( dimensionIndex, axis ); }
	}
	
	if ( bIsCorpusOpen )
	{
		CameraSwitcher ( );
		// TODO - if axis != COLOR, retrain point picker // is this still needed here? already retraining in liveview
	}
}

int ExplorerMenu::GetDimensionIndex ( std::string& dimension )
{
	for ( int i = 0; i < mRawView->GetDimensions ( ).size ( ); i++ )
	{
		if ( mRawView->GetDimensions ( )[i] == dimension )
		{
			return i;
		}
	}
	ofLogWarning ( "LiveView" ) << "Dimension " << dimension << " not found";
	return -1;
}

void ExplorerMenu::CameraSwitcher ( )
{
	bool isXNone = mDimensionDropdownX->getAllSelected ( )[0] == "None";
	bool isYNone = mDimensionDropdownY->getAllSelected ( )[0] == "None";
	bool isZNone = mDimensionDropdownZ->getAllSelected ( )[0] == "None";
	int numDisabledAxes = isXNone + isYNone + isZNone;

	Utils::Axis							  disabledAxis = Utils::Axis::NONE;
	if		( isXNone )					{ disabledAxis = Utils::Axis::X; }
	else if ( isYNone )					{ disabledAxis = Utils::Axis::Y; }
	else if ( isZNone )					{ disabledAxis = Utils::Axis::Z; }
	else if ( numDisabledAxes > 1 )		{ disabledAxis = Utils::Axis::MULTIPLE; }

	bool current3D = mLiveView.Is3D ( );

	if ( disabledAxis == Utils::Axis::NONE || disabledAxis == Utils::Axis::MULTIPLE )
	{
		if ( !mLiveView.Is3D ( ) )
		{
			mLiveView.Set3D ( true );
			mLiveView.Init3DCam ( );
		}
	}
	else
	{
		if ( mLiveView.Is3D ( ) || disabledAxis != mDisabledAxis )
		{
			mLiveView.Set3D ( false );
			mLiveView.Init2DCam ( disabledAxis );
			mDisabledAxis = disabledAxis;
		}
	}
}

// Listener Functions --------------------------

void ExplorerMenu::WindowResized ( )
{
	mMainPanel.setPosition ( ofGetWidth ( ) - mLayout.explorePanelWidth, mLayout.explorePanelOriginY );
	mMorphPanel.setPosition(ofGetWidth() - mLayout.explorePanelWidth, 
		mMainPanel.getPosition().y + mMainPanel.getHeight() + 20);
	mCrossoverPanel.setPosition(ofGetWidth() - mLayout.explorePanelWidth,
		mMorphPanel.getPosition().y + mMorphPanel.getHeight() + 20);
	mAudioDevicePanel.setPosition(ofGetWidth() - mLayout.explorePanelWidth,
		mCrossoverPanel.getPosition().y + mCrossoverPanel.getHeight() + 20);
}

void ExplorerMenu::SwapDimensionX ( string& dimension )
{
	SwapDimension ( dimension, Utils::Axis::X );
}

void ExplorerMenu::SwapDimensionY ( string& dimension )
{
	SwapDimension ( dimension, Utils::Axis::Y );
}

void ExplorerMenu::SwapDimensionZ ( string& dimension )
{
	SwapDimension ( dimension, Utils::Axis::Z );
}

void ExplorerMenu::SwapDimensionColor ( string& dimension )
{
	SwapDimension ( dimension, Utils::Axis::COLOR );
}

void ExplorerMenu::SwitchColorSpectrum ( bool& fullSpectrum )
{
	if ( fullSpectrum ) { mColorSpectrumSwitcher.setName ( "Color Spectrum: Full" ); }
	else { mColorSpectrumSwitcher.setName ( "Color Spectrum: Red<->Blue" ); }
	mLiveView.SetColorFullSpectrum ( fullSpectrum );
	SwapDimension ( mDimensionDropdownColor->getAllSelected ( )[0], Utils::Axis::COLOR );
}

void ExplorerMenu::getAudioBlock(float* output, int bufferSize, int nChannels)
{
	// Forward the audio request to the central AudioPlayback engine so that TIME-mode audio is heard.
	auto playback = mLiveView.GetAudioPlayback();
	if ( playback )
	{
		// Create a temporary buffer for the engine to fill
		ofSoundBuffer temp;
		temp.allocate( bufferSize, nChannels );
		temp.setSampleRate( playback->GetSampleRate() );
		
		// Ask the playback engine to render into the temp buffer
		playback->audioOut( temp );
		
		// Copy the rendered samples into the output buffer supplied by ofApp
		memcpy( output, temp.getBuffer().data(), bufferSize * nChannels * sizeof( float ) );
	}
	else
	{
		// No playback engine available – output silence
		memset( output, 0, bufferSize * nChannels * sizeof( float ) );
	}
}

// Morph Mode Listener Implementations ---------

void ExplorerMenu::onMorphModeChanged ( bool& enabled )
{
	bMorphMode = enabled;
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetMorphMode ( enabled );
	}
}

void ExplorerMenu::onContinuousMorphModeChanged ( bool& enabled )
{
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		// Set the playback mode based on the toggle state
		if ( enabled )
		{
			mLiveView.GetAudioPlayback ( )->SetPlaybackMode ( Explorer::AudioPlayback::PlaybackMode::CONTINUOUS );
		}
		else
		{
			mLiveView.GetAudioPlayback ( )->SetPlaybackMode ( Explorer::AudioPlayback::PlaybackMode::DISCRETE );
		}
	}
}

void ExplorerMenu::onKDTreeNeighborsChanged ( int& count )
{
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetKDTreeNeighbors ( count );
		ofLogNotice("ExplorerMenu") << "KD-tree neighbors changed to: " << count;
	}
}

void ExplorerMenu::onTransitionDurationChanged ( float& duration )
{
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetMorphTransitionDuration ( duration );
	}
}

void ExplorerMenu::onSTFTSizeChanged ( int& size )
{
	// Ensure STFT size is power of 2
	int validSize = size;
	if (size <= 512) validSize = 512;
	else if (size <= 1024) validSize = 1024;
	else if (size <= 2048) validSize = 2048;
	else validSize = 4096;
	
	if (validSize != size) {
		mSTFTSizeSlider = validSize;
	}
	
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetMorphSTFTSize ( validSize );
	}
}

void ExplorerMenu::onInterpolationPositionChanged ( float& position )
{
	// TODO: This would require manual morph control implementation
	// Currently morphing is automatic based on crossfade progress
}

void ExplorerMenu::onLoopModeChanged ( bool& enabled )
{
	// TODO: This would require loop control for morphing
	// Currently handled by existing loop playhead setting
}

// Crossover Jump Listener Implementations -----

void ExplorerMenu::onLoopPlayheadsChanged ( bool& enabled )
{
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetLoopPlayheads ( enabled );
	}
}

void ExplorerMenu::onJumpSameFileAllowedChanged ( bool& allowed )
{
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetJumpSameFileAllowed ( allowed );
	}
}

void ExplorerMenu::onJumpSameFileMinTimeDiffChanged ( int& timeDiff )
{
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetJumpSameFileMinTimeDiff ( timeDiff );
	}
}

void ExplorerMenu::onCrossoverJumpChanceChanged ( float& chance )
{
	// Convert percentage to jumps per thousand
	int jumpsPerThousand = static_cast<int>(chance * 10.0f);
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetCrossoverJumpChance ( jumpsPerThousand );
	}
}

void ExplorerMenu::onCrossfadeMaxSampleLengthChanged ( int& length )
{
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetCrossfadeSampleLength ( length );
	}
}

void ExplorerMenu::onMaxJumpDistanceSpaceChanged ( float& distance )
{
	// Convert percentage to distance x1000
	int distanceX1000 = static_cast<int>(distance * 10.0f);
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetMaxJumpDistanceSpace ( distanceX1000 );
	}
}

void ExplorerMenu::onMaxJumpTargetsChanged ( int& targets )
{
	if ( mLiveView.GetAudioPlayback ( ) )
	{
		mLiveView.GetAudioPlayback ( )->SetMaxJumpTargets ( targets );
	}
}

void ExplorerMenu::onBufferSizeChanged ( string& value )
{
	// Convert string to int and find index
	int bufferSize = std::stoi(value);
	int index = -1;
	
	for (int i = 0; i < mBufferSizes.size(); i++) {
		if (mBufferSizes[i] == bufferSize) {
			index = i;
			break;
		}
	}
	
	if ( index >= 0 && index < mBufferSizes.size() && mLiveView.GetAudioPlayback ( ) )
	{
		mCurrentBufferSizeIndex = index;
		
		// Get current sample rate and device
		int sampleRate = 44100; // Default, should get from current stream
		ofSoundDevice device = mCurrentDeviceIndex < mAudioDevices.size() ? 
			mAudioDevices[mCurrentDeviceIndex] : ofSoundDevice();
		
		// Restart audio with new buffer size
		mLiveView.GetAudioPlayback ( )->RestartAudio ( sampleRate, bufferSize, device );
	}
}

void ExplorerMenu::onOutDeviceChanged ( string& deviceName )
{
	// Find the device with matching name
	int actualDeviceIndex = -1;
	
	for ( int i = 0; i < mAudioDevices.size(); i++ )
	{
		if ( mAudioDevices[i].name == deviceName && mAudioDevices[i].outputChannels > 0 )
		{
			actualDeviceIndex = i;
			break;
		}
	}
	
	if ( actualDeviceIndex >= 0 && actualDeviceIndex < mAudioDevices.size() && mLiveView.GetAudioPlayback ( ) )
	{
		mCurrentDeviceIndex = actualDeviceIndex;
		
		// Get current buffer size and sample rate
		int bufferSize = mBufferSizes[mCurrentBufferSizeIndex];
		int sampleRate = 44100; // Default, should get from current stream
		
		// Restart audio with new device
		mLiveView.GetAudioPlayback ( )->RestartAudio ( sampleRate, bufferSize, mAudioDevices[actualDeviceIndex] );
	}
}
