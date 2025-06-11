#include "./PointPicker.h"
#include "./SpaceDefs.h"

#include <ofGraphics.h>
#include <of3DGraphics.h>
#include <queue>
#include <algorithm>
#include <cmath>

using namespace Acorex;

void Explorer::PointPicker::Initialise ( const Utils::DataSet& dataset, const Utils::DimensionBounds& dimensionBounds )
{
	mFullFluidSet = fluid::FluidDataSet<std::string, double, 1> ( dataset.dimensionNames.size ( ) );
	mLiveFluidSet = fluid::FluidDataSet<std::string, double, 1> ( 3 );


	bPicker = false; b3D = true; bSkipTraining = true; bNearestCheckNeeded = false;

	bDimensionsFilled[0] = false; bDimensionsFilled[1] = false; bDimensionsFilled[2] = false;
	mDimensionsIndices[0] = -1; mDimensionsIndices[1] = -1; mDimensionsIndices[2] = -1;

	mNearestPoint = -1; mNearestDistance = -1; mNearestPointFile = -1; mNearestPointTime = -1;

	if ( mCorpusFileLookUp.size ( ) > 0 ) { mCorpusFileLookUp.clear ( ); }
	if ( mCorpusTimeLookUp.size ( ) > 0 ) { mCorpusTimeLookUp.clear ( ); }

	Utils::DataSet scaledDataset = dataset;
	ScaleDataset ( scaledDataset, dimensionBounds );

	for ( int file = 0; file < dataset.fileList.size ( ); file++ )
	{
		if ( dataset.analysisSettings.bTime )
		{
			for ( int timepoint = 0; timepoint < dataset.time.raw[file].size ( ); timepoint++ )
			{
				mCorpusFileLookUp.push_back ( file );
				mCorpusTimeLookUp.push_back ( timepoint );
			}
		}
		else
		{
			mCorpusFileLookUp.push_back ( file );
			mCorpusTimeLookUp.push_back ( 0 ); // In STATS mode, use time index 0
		}
	}

    std::vector<int> tempVector;
    
	mDatasetConversion.CorpusToFluid ( mFullFluidSet, scaledDataset, tempVector );

	if ( !bListenersAdded )
	{
		ofAddListener ( ofEvents ( ).mouseMoved, this, &Explorer::PointPicker::MouseMoved );
		ofAddListener ( ofEvents ( ).keyReleased, this, &Explorer::PointPicker::KeyEvent );
		ofAddListener ( ofEvents ( ).mouseReleased, this, &Explorer::PointPicker::MouseReleased );
		bListenersAdded = true;
	}
}

void Explorer::PointPicker::Train ( int dimensionIndex, Utils::Axis axis, bool none )
{
	if ( axis == Utils::Axis::X ) { bDimensionsFilled[0] = !none; mDimensionsIndices[0] = dimensionIndex; }
	else if ( axis == Utils::Axis::Y ) { bDimensionsFilled[1] = !none; mDimensionsIndices[1] = dimensionIndex; }
	else if ( axis == Utils::Axis::Z ) { bDimensionsFilled[2] = !none; mDimensionsIndices[2] = dimensionIndex; }
	else { return; }

	int dimsFilled = bDimensionsFilled[0] + bDimensionsFilled[1] + bDimensionsFilled[2];
	if ( dimsFilled < 2 ) { bTrained = false; return; }

	if ( axis == Utils::Axis::Z ) { bSkipTraining = false; }
	if ( bSkipTraining ) { return; }

	mLiveFluidSet = fluid::FluidDataSet<std::string, double, 1> ( dimsFilled );

	for ( int point = 0; point < mFullFluidSet.size ( ); point++ )
	{
		fluid::RealVector pointData ( dimsFilled );
		if ( dimsFilled == 3 || bDimensionsFilled[2] == false )
		{
			for ( int dim = 0; dim < dimsFilled; dim++ )
			{
				pointData[dim] = mFullFluidSet.get ( mFullFluidSet.getIds ( )[point] )[mDimensionsIndices[dim]];
			}
		}
		else if ( bDimensionsFilled[1] == false )
		{
			pointData[0] = mFullFluidSet.get ( mFullFluidSet.getIds ( )[point] )[mDimensionsIndices[0]];
			pointData[1] = mFullFluidSet.get ( mFullFluidSet.getIds ( )[point] )[mDimensionsIndices[2]];
		}
		else if ( bDimensionsFilled[0] == false )
		{
			pointData[0] = mFullFluidSet.get ( mFullFluidSet.getIds ( )[point] )[mDimensionsIndices[1]];
			pointData[1] = mFullFluidSet.get ( mFullFluidSet.getIds ( )[point] )[mDimensionsIndices[2]];
		}

		mLiveFluidSet.add ( mFullFluidSet.getIds ( )[point], pointData );
	}

	ofLogNotice ( "PointPicker" ) << "Training KDTree...";
	mKDTree = fluid::algorithm::KDTree ( mLiveFluidSet );
	ofLogNotice ( "PointPicker" ) << "KDTree Trained.";
	bTrained = true;

	if ( dimsFilled == 2 ) { b3D = false; }
	if ( dimsFilled == 3 ) { b3D = true; }
}

void Explorer::PointPicker::Exit ( )
{
	RemoveListeners ( );
}

void Explorer::PointPicker::RemoveListeners ( )
{
	if ( bListenersAdded )
	{
		ofRemoveListener ( ofEvents ( ).mouseMoved, this, &Explorer::PointPicker::MouseMoved );
		ofRemoveListener ( ofEvents ( ).keyReleased, this, &Explorer::PointPicker::KeyEvent );
		ofRemoveListener ( ofEvents ( ).mouseReleased, this, &Explorer::PointPicker::MouseReleased );
		bListenersAdded = false;
	}
}

void Explorer::PointPicker::ScaleDataset ( Utils::DataSet& scaledDataset, const Utils::DimensionBounds& dimensionBounds )
{
	if ( scaledDataset.analysisSettings.bTime )
	{
		for ( int file = 0; file < scaledDataset.time.raw.size ( ); file++ )
		{
			for ( int timepoint = 0; timepoint < scaledDataset.time.raw[file].size ( ); timepoint++ )
			{
				for ( int dimension = 0; dimension < scaledDataset.dimensionNames.size ( ); dimension++ )
				{
					scaledDataset.time.raw[file][timepoint][dimension] = ofMap ( 
						scaledDataset.time.raw[file][timepoint][dimension],
						dimensionBounds.GetMinBound ( dimension ),
						dimensionBounds.GetMaxBound ( dimension ),
						0.0, 1.0, false );
				}
			}
		}
	}
	else if ( scaledDataset.analysisSettings.hasBeenReduced )
	{
		for ( int file = 0; file < scaledDataset.stats.reduced.size ( ); file++ )
		{
			for ( int dimension = 0; dimension < scaledDataset.dimensionNames.size ( ); dimension++ )
			{
				scaledDataset.stats.reduced[file][dimension] = ofMap (
					scaledDataset.stats.reduced[file][dimension],
					dimensionBounds.GetMinBound ( dimension ),
					dimensionBounds.GetMaxBound ( dimension ),
					0.0, 1.0, false );
			}
		}
	}
	else
	{
		for ( int file = 0; file < scaledDataset.stats.raw.size ( ); file++ )
		{
			for ( int dimension = 0; dimension < scaledDataset.stats.raw[file].size ( ); dimension++ )
			{
				int statistic = dimension % DATA_NUM_STATS;
				int dividedDimension = dimension / DATA_NUM_STATS;
				scaledDataset.stats.raw[file][dividedDimension][statistic] = ofMap (
					scaledDataset.stats.raw[file][dividedDimension][statistic],
					dimensionBounds.GetMinBound ( dimension ),
					dimensionBounds.GetMaxBound ( dimension ),
					0.0, 1.0, false );
			}
		}
	}
}

void Explorer::PointPicker::SlowUpdate ( )
{
	FindNearest ( );
}

void Explorer::PointPicker::Draw ( )
{
	if ( bDebug )
	{
		if ( mNearestPoint != -1 )
		{
			ofDrawBitmapStringHighlight ( "Nearest Point: " + std::to_string ( mNearestPoint ), 20, ofGetHeight ( ) - 100 );
			ofDrawBitmapStringHighlight ( "Nearest Distance: " + std::to_string ( mNearestDistance ), 20, ofGetHeight ( ) - 80 );
		}

		ofEnableDepthTest ( );
		mCamera->begin ( );

		ofSetColor ( 150, 150, 255, 125 );
		for ( int i = 0; i < testPoints.size ( ); i++ )
		{
			ofDrawSphere ( testPoints[i], testRadii[i] );
		}

		ofSetColor ( 255, 255, 255, 25 );
		for ( int i = 0; i < testPointsOutOfRange.size ( ); i++ )
		{
			ofDrawSphere ( testPointsOutOfRange[i], testRadiiOutOfRange[i] );
		}

		mCamera->end ( );
		ofDisableDepthTest ( );
	}
}

void Explorer::PointPicker::FindNearest ( )
{
	if ( !ofGetMousePressed ( 2 ) && !bPicker ) { return; }
	if ( !bTrained ) { return; }
	if ( !bNearestCheckNeeded ) { return; }
	bNearestCheckNeeded = false;

	mNearestPoint = -1; mNearestPointFile = -1; mNearestPointTime = -1;
	mNearestDistance = std::numeric_limits<double>::max ( );

	int mouseX = ofGetMouseX ( );
	int mouseY = ofGetMouseY ( );

	if ( !b3D )
	{
		// 2D nearest

		glm::vec3 rayPosition = mCamera->screenToWorld ( glm::vec3 ( mouseX, mouseY, 0 ) );
		glm::vec2 rayPosition2D;
		
		if ( !bDimensionsFilled[0] ) { rayPosition2D.x = rayPosition.y; rayPosition2D.y = rayPosition.z; rayPosition.x = 0; }
		if ( !bDimensionsFilled[1] ) { rayPosition2D.x = rayPosition.x; rayPosition2D.y = rayPosition.z; rayPosition.y = 0; }
		if ( !bDimensionsFilled[2] ) { rayPosition2D.x = rayPosition.x; rayPosition2D.y = rayPosition.y; rayPosition.z = 0; }

		fluid::RealVector query ( 2 );

		query[0] = ofMap ( rayPosition2D.x, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false );
		query[1] = ofMap ( rayPosition2D.y, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false );

		double maxAllowedDistance = ofMap ( mCamera->getScale ( ).x, SpaceDefs::mZoomMin2D, SpaceDefs::mZoomMax2D, maxAllowedDistanceNear * 1.5, maxAllowedDistanceFar * 1.5 );

		auto [dist, id] = mKDTree.kNearest ( query, 1, maxAllowedDistance );

		if ( dist.size ( ) == 0 ) { return; }

		if ( dist[0] < mNearestDistance )
		{
			mNearestDistance = dist[0];
			mNearestPoint = std::stoi ( *id[0] );
			mNearestPointFile = mCorpusFileLookUp[mNearestPoint];
			if ( mCorpusTimeLookUp.size ( ) > 0 ) { mNearestPointTime = mCorpusTimeLookUp[mNearestPoint]; }
		}

		return;
	}

	// 3D nearest
	
	double desiredRayLength = 15000.0f;
	double rayLength = 0.0f;
	std::vector<double> rayPointSpacing;
	
	do
	{
		double maxAllowedDistance = ofMap ( rayLength, 0.0f, desiredRayLength, 0.01, 0.05, false );
		rayPointSpacing.push_back ( maxAllowedDistance );
		rayLength += maxAllowedDistance * ( SpaceDefs::mSpaceMax - SpaceDefs::mSpaceMin );
	} while ( rayLength < desiredRayLength );

	glm::vec3 rayDirection = mCamera->screenToWorld ( glm::vec3 ( (float)mouseX, (float)mouseY, 0.0f ) );
	rayDirection = glm::normalize ( rayDirection - mCamera->getPosition ( ) );
	double depth = 0.0f;

	if ( testPoints.size ( ) > 0 ) { testPoints.clear ( ); }
	if ( testRadii.size ( ) > 0 ) { testRadii.clear ( ); }
	if ( testPointsOutOfRange.size ( ) > 0 ) { testPointsOutOfRange.clear ( ); }
	if ( testRadiiOutOfRange.size ( ) > 0 ) { testRadiiOutOfRange.clear ( ); }

	for ( int rayPoint = 1; rayPoint < rayPointSpacing.size ( ); rayPoint++ )
	{
		depth += rayPointSpacing[rayPoint] * ( SpaceDefs::mSpaceMax - SpaceDefs::mSpaceMin );

		glm::vec3 rayPointPosition = mCamera->getPosition ( ) + glm::vec3 ( rayDirection.x * depth,
																			rayDirection.y * depth, 
																			rayDirection.z * depth );

		if ( rayPointPosition.x < ( SpaceDefs::mSpaceMin - 250 ) || rayPointPosition.x > ( SpaceDefs::mSpaceMax + 250 ) ||
			rayPointPosition.y < ( SpaceDefs::mSpaceMin - 250 ) || rayPointPosition.y > ( SpaceDefs::mSpaceMax + 250 ) ||
			rayPointPosition.z < ( SpaceDefs::mSpaceMin - 250 ) || rayPointPosition.z > ( SpaceDefs::mSpaceMax + 250 ) )
		{
			if ( bDebug )
			{
				testPointsOutOfRange.push_back ( rayPointPosition );
				testRadiiOutOfRange.push_back ( rayPointSpacing[rayPoint] * ( SpaceDefs::mSpaceMax - SpaceDefs::mSpaceMin ) );
			}
			continue; 
		}

		if ( bDebug )
		{
			testPoints.push_back ( rayPointPosition );
			testRadii.push_back ( rayPointSpacing[rayPoint] * (SpaceDefs::mSpaceMax - SpaceDefs::mSpaceMin) );
		}

		fluid::RealVector query ( 3 );

		query[0] = ofMap ( rayPointPosition.x, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false );
		query[1] = ofMap ( rayPointPosition.y, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false );
		query[2] = ofMap ( rayPointPosition.z, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false );

		auto [dist, id] = mKDTree.kNearest ( query, 1, rayPointSpacing[rayPoint] );

		if ( dist.size ( ) == 0 ) { continue; }

		if ( dist[0] < mNearestDistance )
		{
			mNearestDistance = dist[0];
			mNearestPoint = std::stoi ( *id[0] );
			mNearestPointFile = mCorpusFileLookUp[mNearestPoint];
			if ( mCorpusTimeLookUp.size ( ) > 0 ) { mNearestPointTime = mCorpusTimeLookUp[mNearestPoint]; }
		}
	}
}

void Explorer::PointPicker::KeyEvent ( ofKeyEventArgs& args )
{
	if ( args.type == ofKeyEventArgs::Type::Released )
	{
		if ( args.key == OF_KEY_F3 ) { bDebug = !bDebug; }
		else if ( args.key == OF_KEY_TAB ) { bPicker = !bPicker; }
	}
}


void Explorer::PointPicker::FindNearestToMouse ( )
{
	FindNearest();
}

bool Explorer::PointPicker::FindNearestToPosition ( const glm::vec3& position, Utils::PointFT& nearestPoint, Utils::PointFT currentPoint, 
                                int maxAllowedDistanceSpaceX1000, int maxAllowedTargets, bool sameFileAllowed, 
                                int minTimeDiffSameFile, int remainingSamplesRequired, const Utils::AudioData& audioSet, size_t hopSize )
{
    if ( maxAllowedDistanceSpaceX1000 == 0 ) { return false; }

    if ( mPointPickerMutex.try_lock ( ) )
    {
        std::lock_guard<std::mutex> lock ( mPointPickerMutex, std::adopt_lock );

        double maxAllowedDistanceSpace = static_cast<double>( maxAllowedDistanceSpaceX1000 ) / 1000.0;

        if ( !b3D )
        {
            // ---------------- 2D MODE ----------------
            glm::vec2 position2D;

            if ( !bDimensionsFilled[0] ) { position2D.x = position.y; position2D.y = position.z; }
            if ( !bDimensionsFilled[1] ) { position2D.x = position.x; position2D.y = position.z; }
            if ( !bDimensionsFilled[2] ) { position2D.x = position.x; position2D.y = position.y; }

            fluid::RealVector query ( 2 );
            query[0] = ofMap ( position2D.x, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false );
            query[1] = ofMap ( position2D.y, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false );

            auto [dist, id] = mKDTree.kNearest ( query, maxAllowedTargets, maxAllowedDistanceSpace );
            if ( dist.size ( ) == 0 ) { return false; }

            double nearestDistance = std::numeric_limits<double>::max ( );
            bool   jumpFound       = false;

            for ( size_t i = 0; i < dist.size ( ); ++i )
            {
                if ( dist[i] >= nearestDistance ) { continue; }

                int point = std::stoi ( *id[i] );

                if ( !sameFileAllowed && mCorpusFileLookUp[point] == currentPoint.file ) { continue; }
                size_t timeDiff = mCorpusTimeLookUp[point] > currentPoint.time ?
                                   mCorpusTimeLookUp[point] - currentPoint.time :
                                   currentPoint.time - mCorpusTimeLookUp[point];
                if ( sameFileAllowed && mCorpusFileLookUp[point] == currentPoint.file && timeDiff < minTimeDiffSameFile ) { continue; }

                if ( audioSet.raw[mCorpusFileLookUp[point]].getNumFrames ( ) - ( (size_t)mCorpusTimeLookUp[point] * hopSize ) < remainingSamplesRequired )
                {
                    continue; // not enough samples left
                }

                nearestDistance = dist[i];
                nearestPoint.file = mCorpusFileLookUp[point];
                nearestPoint.time = mCorpusTimeLookUp[point];
                jumpFound = true;
            }

            return jumpFound;
        }

        // ---------------- 3D MODE ----------------
        fluid::RealVector query ( 3 );
        query[0] = ofMap ( position.x, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false );
        query[1] = ofMap ( position.y, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false );
        query[2] = ofMap ( position.z, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false );

        auto [dist, id] = mKDTree.kNearest ( query, maxAllowedTargets, maxAllowedDistanceSpace );
        if ( dist.size ( ) == 0 ) { return false; }

        double nearestDistance = std::numeric_limits<double>::max ( );
        bool   jumpFound       = false;

        for ( size_t i = 0; i < dist.size ( ); ++i )
        {
            if ( dist[i] >= nearestDistance ) { continue; }

            int point = std::stoi ( *id[i] );

            if ( !sameFileAllowed && mCorpusFileLookUp[point] == currentPoint.file ) { continue; }
            size_t timeDiff = mCorpusTimeLookUp[point] > currentPoint.time ?
                               mCorpusTimeLookUp[point] - currentPoint.time :
                               currentPoint.time - mCorpusTimeLookUp[point];
            if ( sameFileAllowed && mCorpusFileLookUp[point] == currentPoint.file && timeDiff < minTimeDiffSameFile ) { continue; }

            if ( audioSet.raw[mCorpusFileLookUp[point]].getNumFrames ( ) - ( (size_t)mCorpusTimeLookUp[point] * hopSize ) < remainingSamplesRequired )
            {
                continue; // not enough samples left
            }

            nearestDistance = dist[i];
            nearestPoint.file = mCorpusFileLookUp[point];
            nearestPoint.time = mCorpusTimeLookUp[point];
            jumpFound = true;
        }

        return jumpFound;
    }

    return false;
}

void Explorer::PointPicker::MouseReleased ( ofMouseEventArgs& args )
{
    if ( args.button == 2 ) { bClicked = true; }
}

bool Explorer::PointPicker::FindKNearestToPosition(const glm::vec3& position, KNNResult& result, int k, 
                                                   double maxAllowedDistanceSpace, Utils::PointFT currentPoint, 
                                                   bool sameFileAllowed, int minTimeDiffSameFile, 
                                                   int remainingSamplesRequired, const Utils::AudioData& audioSet, size_t hopSize)
{
    result.clear();
    
    if (!bTrained || k <= 0) { return false; }
    
    if (mPointPickerMutex.try_lock())
    {
        std::lock_guard<std::mutex> lock(mPointPickerMutex, std::adopt_lock);
        
        fluid::RealVector query;
        
        if (!b3D)
        {
            // 2D mode
            glm::vec2 position2D;
            
            if (!bDimensionsFilled[0]) { position2D.x = position.y; position2D.y = position.z; }
            if (!bDimensionsFilled[1]) { position2D.x = position.x; position2D.y = position.z; }
            if (!bDimensionsFilled[2]) { position2D.x = position.x; position2D.y = position.y; }
            
            query = fluid::RealVector(2);
            query[0] = ofMap(position2D.x, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false);
            query[1] = ofMap(position2D.y, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false);
        }
        else
        {
            // 3D mode
            query = fluid::RealVector(3);
            query[0] = ofMap(position.x, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false);
            query[1] = ofMap(position.y, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false);
            query[2] = ofMap(position.z, SpaceDefs::mSpaceMin, SpaceDefs::mSpaceMax, 0.0, 1.0, false);
        }
        
        // Get more candidates than k to account for filtering
        int candidateK = std::min(k * 3, (int)mKDTree.size());
        auto [dist, id] = mKDTree.kNearest(query, candidateK, maxAllowedDistanceSpace);
        
        if (dist.size() == 0) { return false; }
        
        // Process candidates and filter based on constraints
        for (size_t i = 0; i < dist.size() && result.size() < k; ++i)
        {
            int point = std::stoi(*id[i]);
            int fileIdx = mCorpusFileLookUp[point];
            int timeIdx = mCorpusTimeLookUp[point];
            
            // Check file constraint
            if (!sameFileAllowed && fileIdx == currentPoint.file) { continue; }
            
            // Check time difference constraint
            if (sameFileAllowed && fileIdx == currentPoint.file)
            {
                size_t timeDiff = timeIdx > currentPoint.time ?
                                  timeIdx - currentPoint.time :
                                  currentPoint.time - timeIdx;
                if (timeDiff < minTimeDiffSameFile) { continue; }
            }
            
            // Check remaining samples constraint
            if (audioSet.raw[fileIdx].getNumFrames() - ((size_t)timeIdx * hopSize) < remainingSamplesRequired)
            {
                continue;
            }
            
            // Add valid point to result
            Utils::PointFT pt;
            pt.file = fileIdx;
            pt.time = timeIdx;
            
            result.points.push_back(pt);
            result.distances.push_back(dist[i]);
        }
        
        if (result.size() == 0) { return false; }
        
        // Calculate normalized weights using inverse distance weighting
        result.weights.resize(result.size());
        double totalWeight = 0.0;
        
        // Handle edge case where a point has zero distance
        bool hasZeroDistance = false;
        size_t zeroDistIndex = 0;
        
        for (size_t i = 0; i < result.distances.size(); ++i)
        {
            if (result.distances[i] < 1e-10) // Effectively zero
            {
                hasZeroDistance = true;
                zeroDistIndex = i;
                break;
            }
        }
        
        if (hasZeroDistance)
        {
            // If we have a zero distance point, give it all the weight
            std::fill(result.weights.begin(), result.weights.end(), 0.0);
            result.weights[zeroDistIndex] = 1.0;
        }
        else
        {
            // Calculate inverse distance weights
            for (size_t i = 0; i < result.distances.size(); ++i)
            {
                result.weights[i] = 1.0 / result.distances[i];
                totalWeight += result.weights[i];
            }
            
            // Normalize weights to sum to 1.0
            if (totalWeight > 0.0)
            {
                for (size_t i = 0; i < result.weights.size(); ++i)
                {
                    result.weights[i] /= totalWeight;
                }
            }
        }
        
        return true;
    }
    
    return false;
}

void Explorer::PointPicker::FindKNearestToMouse(KNNResult& result, int k)
{
    if (!bTrained) { return; }
    
    int mouseX = ofGetMouseX();
    int mouseY = ofGetMouseY();
    
    glm::vec3 position;
    
    if (!b3D)
    {
        // 2D mode - use screen to world conversion
        position = mCamera->screenToWorld(glm::vec3(mouseX, mouseY, 0));
    }
    else
    {
        // 3D mode - use ray casting to find position
        glm::vec3 rayDirection = mCamera->screenToWorld(glm::vec3((float)mouseX, (float)mouseY, 0.0f));
        rayDirection = glm::normalize(rayDirection - mCamera->getPosition());
        
        // Use a default depth for 3D picking
        float depth = 1000.0f;
        position = mCamera->getPosition() + rayDirection * depth;
    }
    
    // Use a reasonable max distance based on zoom
    double maxDist = b3D ? 0.1 : ofMap(mCamera->getScale().x, SpaceDefs::mZoomMin2D, SpaceDefs::mZoomMax2D, 
                                        maxAllowedDistanceNear * 1.5, maxAllowedDistanceFar * 1.5);
    
    // Create a dummy current point (not used for mouse picking)
    Utils::PointFT dummyPoint;
    dummyPoint.file = -1;
    dummyPoint.time = -1;
    
    // Create a dummy audio set (not used for mouse picking)
    Utils::AudioData dummyAudioSet;
    
    FindKNearestToPosition(position, result, k, maxDist, dummyPoint, true, 0, 0, dummyAudioSet, 1);
}
