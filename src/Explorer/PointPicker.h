#pragma once

#include "./RawView.h"
#include "Utils/DimensionBounds.h"
#include "Utils/DatasetConversion.h"
#include <algorithms/public/KDTree.hpp>
#include <data/FluidDataSet.hpp>
#include <ofCamera.h>
#include <ofEvents.h>
#include <mutex>
#include <vector>

namespace Acorex {
namespace Explorer {

struct KNNResult {
	std::vector<Utils::PointFT> points;
	std::vector<double> distances;
	std::vector<double> weights;
	
	void clear() {
		points.clear();
		distances.clear();
		weights.clear();
	}
	
	size_t size() const { return points.size(); }
};

class PointPicker {
public:
	PointPicker ( ) { }
	~PointPicker ( ) { }

	void Initialise ( const Utils::DataSet& dataset, const Utils::DimensionBounds& dimensionBounds );

	void Train ( int dimensionIndex, Utils::Axis axis, bool none );

	void Exit ( );
	void RemoveListeners ( );

	void SlowUpdate ( );
	void Draw ( );

	void FindNearest ( );
	void FindNearestToMouse ( );
	bool FindNearestToPosition (const glm::vec3& position, Utils::PointFT& nearestPoint, Utils::PointFT currentPoint, 
								int maxAllowedDistanceSpaceX1000, int maxAllowedTargets, bool sameFileAllowed, 
								int minTimeDiffSameFile, int remainingSamplesRequired, const Utils::AudioData& audioSet, size_t hopSize );
	
	// k-NN methods
	bool FindKNearestToPosition(const glm::vec3& position, KNNResult& result, int k, 
								double maxAllowedDistanceSpace, Utils::PointFT currentPoint, 
								bool sameFileAllowed, int minTimeDiffSameFile, 
								int remainingSamplesRequired, const Utils::AudioData& audioSet, size_t hopSize);
	
	void FindKNearestToMouse(KNNResult& result, int k = 3);

	// Setters & Getters ----------------------------

	void SetCamera ( std::shared_ptr<ofCamera> camera ) { mCamera = camera; }
	void SetNearestCheckNeeded ( ) { bNearestCheckNeeded = true; }

	int GetNearestMousePointFile ( ) const { return mNearestPointFile; }
	int GetNearestMousePointTime ( ) const { return mNearestPointTime; }
	double GetNearestMouseDistance ( ) const { return mNearestDistance; }
	int GetNearestPointFile ( ) const { return mNearestPointFile; }
	int GetNearestPointTime ( ) const { return mNearestPointTime; }
	double GetNearestDistance ( ) const { return mNearestDistance; }
	bool IsTrained ( ) const { return bTrained; }
	std::vector<int>* GetCorpusFileLookUp ( ) { return &mCorpusFileLookUp; }
	std::vector<int>* GetCorpusTimeLookUp ( ) { return &mCorpusTimeLookUp; }

private:
	void ScaleDataset ( Utils::DataSet& scaledDataset, const Utils::DimensionBounds& dimensionBounds );

	// Listeners ------------------------------------

	void MouseMoved ( ofMouseEventArgs& args ) { bNearestCheckNeeded = true; }
	void KeyEvent ( ofKeyEventArgs& args );
	void MouseReleased ( ofMouseEventArgs& args );

	// States ---------------------------------------

	bool bPicker = false;
	bool bDebug = false;
	bool bClicked = false;
	bool bDraw = false;
	bool b3D = true;
	bool bTrained = false;
	bool bSkipTraining = true;
	bool bListenersAdded = false;
	bool bNearestCheckNeeded = false;
	bool bDimensionsFilled[3] = { false, false, false };

	// Variables ------------------------------------

	std::shared_ptr<ofCamera> mCamera;

	int mDimensionsIndices[3] = { -1, -1, -1 };

	int mNearestPoint = -1;
	double mNearestDistance = -1;

	double maxAllowedDistanceFar = 0.05;
	double maxAllowedDistanceNear = 0.01;
	
	// k-NN configuration
	static constexpr int DEFAULT_K = 3;

	fluid::algorithm::KDTree mKDTree;

	fluid::FluidDataSet<std::string, double, 1> mFullFluidSet;
	fluid::FluidDataSet<std::string, double, 1> mLiveFluidSet;
	std::vector<int> mCorpusFileLookUp; int mNearestPointFile = -1;
	std::vector<int> mCorpusTimeLookUp; int mNearestPointTime = -1;
	
	Utils::DatasetConversion mDatasetConversion;

	std::vector<glm::vec3> testPoints;
	std::vector<float> testRadii;
	std::vector<glm::vec3> testPointsOutOfRange;
	std::vector<float> testRadiiOutOfRange;

	// Thread safety --------------------------------
	std::mutex mPointPickerMutex;
};

} // namespace Explorer
} // namespace Acorex