#pragma once

#include "Explorer/RawView.h"
#include "Explorer/LiveView.h"
#include "Utils/InterfaceDefs.h"
#include <ofxGui.h>
#include <ofxDropdown.h>

namespace Acorex {

class ExplorerMenu {
public:
	ExplorerMenu ( ) { }
	~ExplorerMenu ( ) { }

	void Initialise ( bool HiDpi );
	void Show ( );
	void Hide ( );
	void Draw ( );
	void Update ( );
	void Exit ( );

	void WindowResized ( );
	
	void getAudioBlock(float* output, int bufferSize, int nChannels);

private:
	void SlowUpdate ( );
	void RemoveListeners ( );

	// Main Functions ------------------------------

	void OpenCorpus ( );
	void SwapDimension ( string dimension, Utils::Axis axis );
	int GetDimensionIndex ( std::string& dimension );
	void CameraSwitcher ( );

	// Listener Functions --------------------------

	void SwapDimensionX ( string& dimension );
	void SwapDimensionY ( string& dimension );
	void SwapDimensionZ ( string& dimension );
	void SwapDimensionColor ( string& dimension );
	void SwitchColorSpectrum ( bool& fullSpectrum );
	
	// Morph Mode Listener Functions ---------------
	void onMorphModeChanged ( bool& enabled );
	void onContinuousMorphModeChanged ( bool& enabled );
	void onKDTreeNeighborsChanged ( int& count );
	void onTransitionDurationChanged ( float& duration );
	void onSTFTSizeChanged ( int& size );
	void onInterpolationPositionChanged ( float& position );
	void onLoopModeChanged ( bool& enabled );
	
	// Crossover Jump Listener Functions -----------
	void onLoopPlayheadsChanged ( bool& enabled );
	void onJumpSameFileAllowedChanged ( bool& allowed );
	void onJumpSameFileMinTimeDiffChanged ( int& timeDiff );
	void onCrossoverJumpChanceChanged ( float& chance );
	void onCrossfadeMaxSampleLengthChanged ( int& length );
	void onMaxJumpDistanceSpaceChanged ( float& distance );
	void onMaxJumpTargetsChanged ( int& targets );
	void onBufferSizeChanged ( string& value );
	void onOutDeviceChanged ( string& deviceName );

	// States --------------------------------------

	bool bDraw = false;

	bool bIsCorpusOpen = false; bool bBlockDimensionFilling = false;
	bool bOpenCorpusDrawWarning = false;
	bool bInitialiseShouldLoad = false;
	bool bListenersAdded = false;
	
	bool bViewPointerShared = false;
	bool bMorphMode = true;

	Utils::Axis mDisabledAxis = Utils::Axis::NONE;

	// Timing --------------------------------------

	int mLastUpdateTime = 0;
	int mSlowUpdateInterval = 100;

	int mOpenCorpusButtonClickTime = 0;
	int mOpenCorpusButtonTimeout = 3000;

	// Panels --------------------------------------

	ofxPanel mMainPanel;
	ofxLabel mCorpusNameLabel;
	ofxButton mOpenCorpusButton;
	unique_ptr<ofxDropdown> mDimensionDropdownX;
	unique_ptr<ofxDropdown> mDimensionDropdownY;
	unique_ptr<ofxDropdown> mDimensionDropdownZ;
	unique_ptr<ofxDropdown> mDimensionDropdownColor;
	ofxToggle mColorSpectrumSwitcher;
	
	// Morph Mode Controls -------------------------
	ofxPanel mMorphPanel;
	ofxToggle mMorphModeToggle;
	ofxToggle mContinuousMorphModeToggle;
	ofxIntField mKDTreeNeighborsField;
	ofxFloatSlider mTransitionDurationSlider;
	ofxIntSlider mSTFTSizeSlider;
	ofxFloatSlider mInterpolationPositionSlider;
	ofxToggle mLoopModeToggle;
	ofxLabel mLoadingLabel;
	
	// Crossover Jump Controls ---------------------
	ofxPanel mCrossoverPanel;
	ofxToggle mLoopPlayheadsToggle;
	ofxToggle mJumpSameFileAllowedToggle;
	ofxIntSlider mJumpSameFileMinTimeDiffSlider;
	ofxFloatSlider mCrossoverJumpChanceSlider;
	ofxIntSlider mCrossfadeMaxSampleLengthSlider;
	ofxFloatSlider mMaxJumpDistanceSpaceSlider;
	ofxIntSlider mMaxJumpTargetsSlider;
	
	// Audio Device Controls -----------------------
	ofxPanel mAudioDevicePanel;
	unique_ptr<ofxDropdown> mBufferSizeDropdown;
	unique_ptr<ofxDropdown> mOutDeviceDropdown;
	
	// Audio Device Tracking -----------------------
	std::vector<int> mBufferSizes = {64, 128, 256, 512, 1024, 2048, 4096};
	int mCurrentBufferSizeIndex = 2; // Default to 256
	std::vector<ofSoundDevice> mAudioDevices;
	int mCurrentDeviceIndex = 0;

	// Acorex Objects ------------------------------

	std::shared_ptr<Explorer::RawView> mRawView;
	Explorer::LiveView mLiveView;
	Utils::Colors mColors;
	Utils::MenuLayout mLayout;
};

} // namespace Acorex