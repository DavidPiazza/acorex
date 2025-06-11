#include "./JSON.h"
#include "MemoryMappedIO.h"
#include "HDF5IO.h"
#include <ofLog.h>
#include <fstream>
#include <filesystem>

using namespace Acorex;

bool Utils::JSON::Write ( const std::string& outputFile, const DataSet& dataset )
{
	try
	{
		std::ofstream file ( outputFile );
		nlohmann::json j = dataset;

		file << j.dump ( 0 ) << std::endl;
		file.close ( );
	}
	catch ( std::exception& e )
	{
		ofLogError ( "JSON" ) << "failed to write output to " << outputFile << " : " << e.what ( );
		return false;
	}

	return true;
}

bool Utils::JSON::Read ( const std::string& inputFile, DataSet& dataset )
{
	try
	{
		std::ifstream file ( inputFile );
		nlohmann::json j;

		file >> j;
		file.close ( );

		dataset = j.template get<DataSet> ( );
	}
	catch ( std::exception& e )
	{
		ofLogError ( "JSON" ) << "failed to read input " << inputFile << " : " << e.what ( );
		return false;
	}

	return true;
}

bool Utils::JSON::Read ( const std::string& inputFile, AnalysisSettings& settings )
{
	try
	{
		std::ifstream file ( inputFile );
		nlohmann::json j;

		file >> j;
		file.close ( );

		settings = j.template get<AnalysisSettings> ( );
	}
	catch ( std::exception& e )
	{
		ofLogError ( "JSON" ) << "failed to read input " << inputFile << " : " << e.what ( );
		return false;
	}

	return true;
}


#ifndef DATA_CHANGE_CHECK_1
#error "data structure changed, please update json serialization"
#endif

void Utils::to_json ( nlohmann::json& j, const DataSet& a )
{
	j = nlohmann::json {	
		TO_J ( currentPointCount),
		TO_J ( dimensionNames ),
		TO_J ( fileList ),
		TO_J ( time.raw ),
		TO_J ( stats.raw ),
		TO_J ( stats.reduced ),
		TO_J_SETTINGS ( currentDimensionCount ),
		TO_J_SETTINGS ( hasBeenReduced ),
		TO_J_SETTINGS ( bTime ),
		TO_J_SETTINGS ( bPitch ),
		TO_J_SETTINGS ( bLoudness ),
		TO_J_SETTINGS ( bShape ),
		TO_J_SETTINGS ( bMFCC ),
		TO_J_SETTINGS ( bTransport ),
		TO_J_SETTINGS ( windowFFTSize ),
		TO_J_SETTINGS ( hopFraction ),
		TO_J_SETTINGS ( nBands ),
		TO_J_SETTINGS ( nCoefs ),
		TO_J_SETTINGS ( minFreq ),
		TO_J_SETTINGS ( maxFreq ),
		TO_J_SETTINGS ( sampleRate ),
		TO_J_SETTINGS ( transportFFTSize ),
		TO_J_SETTINGS ( transportHopFraction ),
		TO_J_SETTINGS ( transportStorePhase ),
		TO_J_SETTINGS ( transportStoreDerivative ) };
	
	// Only include transport data if bTransport is true
	if (a.analysisSettings.bTransport) {
		j["transport"] = a.transport;
	}
}

void Utils::from_json ( const nlohmann::json& j, DataSet& a )
{
	TO_A ( currentPointCount );
	TO_A ( dimensionNames );
	TO_A ( fileList );
	TO_A ( time.raw );
	TO_A ( stats.raw );
	TO_A ( stats.reduced );
	TO_A_SETTINGS ( currentDimensionCount );
	TO_A_SETTINGS ( hasBeenReduced );
	TO_A_SETTINGS ( bTime );
	TO_A_SETTINGS ( bPitch );
	TO_A_SETTINGS ( bLoudness );
	TO_A_SETTINGS ( bShape );
	TO_A_SETTINGS ( bMFCC );
	// Handle optional bTransport field for backward compatibility
	if (j.contains("bTransport")) {
		TO_A_SETTINGS ( bTransport );
	} else {
		a.analysisSettings.bTransport = false; // Default value for old files
	}
	TO_A_SETTINGS ( windowFFTSize );
	TO_A_SETTINGS ( hopFraction );
	TO_A_SETTINGS ( nBands );
	TO_A_SETTINGS ( nCoefs );
	TO_A_SETTINGS ( minFreq );
	TO_A_SETTINGS ( maxFreq );
	// Handle optional sampleRate field for backward compatibility
	if (j.contains("sampleRate")) {
		TO_A_SETTINGS ( sampleRate );
	} else {
		a.analysisSettings.sampleRate = 44100; // Default value
	}
	// Handle optional transport parameters for backward compatibility
	if (j.contains("transportFFTSize")) {
		TO_A_SETTINGS ( transportFFTSize );
	} else {
		a.analysisSettings.transportFFTSize = 2048; // Default value
	}
	if (j.contains("transportHopFraction")) {
		TO_A_SETTINGS ( transportHopFraction );
	} else {
		a.analysisSettings.transportHopFraction = 4; // Default value
	}
	if (j.contains("transportStorePhase")) {
		TO_A_SETTINGS ( transportStorePhase );
	} else {
		a.analysisSettings.transportStorePhase = true; // Default value
	}
	if (j.contains("transportStoreDerivative")) {
		TO_A_SETTINGS ( transportStoreDerivative );
	} else {
		a.analysisSettings.transportStoreDerivative = true; // Default value
	}
	
	// Load transport data if present and bTransport is true
	if (a.analysisSettings.bTransport && j.contains("transport")) {
		j.at("transport").get_to(a.transport);
	}
}

void Utils::to_json ( nlohmann::json& j, const AnalysisSettings& a )
{
	j = nlohmann::json { 
		TO_J ( currentDimensionCount ),
		TO_J ( hasBeenReduced ),
		TO_J ( bTime ),
		TO_J ( bPitch ),
		TO_J ( bLoudness ),
		TO_J ( bShape ),
		TO_J ( bMFCC ),
		TO_J ( bTransport ),
		TO_J ( windowFFTSize ),
		TO_J ( hopFraction ),
		TO_J ( nBands ),
		TO_J ( nCoefs ),
		TO_J ( minFreq ),
		TO_J ( maxFreq ),
		TO_J ( sampleRate ),
		TO_J ( transportFFTSize ),
		TO_J ( transportHopFraction ),
		TO_J ( transportStorePhase ),
		TO_J ( transportStoreDerivative ) };
}

void Utils::from_json ( const nlohmann::json& j, AnalysisSettings& a )
{ 
	TO_A ( currentDimensionCount );
	TO_A ( hasBeenReduced );
	TO_A ( bTime );
	TO_A ( bPitch );
	TO_A ( bLoudness );
	TO_A ( bShape );
	TO_A ( bMFCC );
	// Handle optional bTransport field for backward compatibility
	if (j.contains("bTransport")) {
		TO_A ( bTransport );
	} else {
		a.bTransport = false; // Default value for old files
	}
	TO_A ( windowFFTSize );
	TO_A ( hopFraction );
	TO_A ( nBands );
	TO_A ( nCoefs );
	TO_A ( minFreq );
	TO_A ( maxFreq );
	// Handle optional sampleRate field for backward compatibility
	if (j.contains("sampleRate")) {
		TO_A ( sampleRate );
	} else {
		a.sampleRate = 44100; // Default value
	}
	// Handle optional transport parameters for backward compatibility
	if (j.contains("transportFFTSize")) {
		TO_A ( transportFFTSize );
	} else {
		a.transportFFTSize = 2048; // Default value
	}
	if (j.contains("transportHopFraction")) {
		TO_A ( transportHopFraction );
	} else {
		a.transportHopFraction = 4; // Default value
	}
	if (j.contains("transportStorePhase")) {
		TO_A ( transportStorePhase );
	} else {
		a.transportStorePhase = true; // Default value
	}
	if (j.contains("transportStoreDerivative")) {
		TO_A ( transportStoreDerivative );
	} else {
		a.transportStoreDerivative = true; // Default value
	}
}

void Utils::to_json ( nlohmann::json& j, const TransportFrame& a )
{
	j = nlohmann::json {
		TO_J ( magnitude ),
		TO_J ( phase ),
		TO_J ( dH )
	};
}

void Utils::from_json ( const nlohmann::json& j, TransportFrame& a )
{
	TO_A ( magnitude );
	TO_A ( phase );
	TO_A ( dH );
}

void Utils::to_json ( nlohmann::json& j, const TransportData& a )
{
	j = nlohmann::json {
		TO_J ( frames )
	};
}

void Utils::from_json ( const nlohmann::json& j, TransportData& a )
{
	TO_A ( frames );
}

// Hybrid I/O not implemented - requires MemoryMappedIO
bool Utils::JSON::WriteHybrid ( const std::string& outputFile, const DataSet& dataset )
{
	try
	{
		// Create a copy of dataset without transport data for JSON
		DataSet datasetCopy = dataset;
		datasetCopy.transport.clear();
		
		// Write main data to JSON
		nlohmann::json j = datasetCopy;
		std::ofstream file ( outputFile );
		file << j.dump ( 4 );
		file.close();
		
		// Write transport data to binary file if present
		if (dataset.analysisSettings.bTransport && dataset.transport.fileCount() > 0) {
			std::filesystem::path jsonPath(outputFile);
			std::filesystem::path binaryPath = jsonPath;
			binaryPath.replace_extension(".transport");
			
			TransportDataIO io;
			if (!io.write(binaryPath.string(), dataset.transport)) {
				// Clean up JSON file on failure
				std::filesystem::remove(outputFile);
				ofLogError("JSON::WriteHybrid") << "Failed to write transport data to: " << binaryPath.string();
				return false;
			}
			
			ofLogNotice("JSON::WriteHybrid") << "Wrote transport data to: " << binaryPath.string();
		}
		
		return true;
	}
	catch ( const std::exception& e )
	{
		ofLogError ( "JSON::WriteHybrid" ) << "Exception writing file: " << e.what ( );
		return false;
	}
}

bool Utils::JSON::ReadHybrid ( const std::string& inputFile, DataSet& dataset )
{
	try
	{
		// Read main data from JSON
		std::ifstream i ( inputFile );
		nlohmann::json j;
		i >> j;
		dataset = j;
		
		// Read transport data from binary file if expected
		if (dataset.analysisSettings.bTransport) {
			std::filesystem::path jsonPath(inputFile);
			std::filesystem::path binaryPath = jsonPath;
			binaryPath.replace_extension(".transport");
			
			if (std::filesystem::exists(binaryPath)) {
				TransportDataIO io;
				if (!io.read(binaryPath.string(), dataset.transport)) {
					// Continue without transport data, but log warning
					ofLogWarning("JSON::ReadHybrid") << "Failed to load transport data from: " << binaryPath.string();
					dataset.transport.clear();
				} else {
					ofLogNotice("JSON::ReadHybrid") << "Loaded transport data from: " << binaryPath.string();
				}
			} else {
				ofLogWarning("JSON::ReadHybrid") << "Transport data file not found: " << binaryPath.string();
			}
		}
		
		return true;
	}
	catch ( const std::exception& e )
	{
		ofLogError ( "JSON::ReadHybrid" ) << "Exception reading file: " << e.what ( );
		return false;
	}
}

bool Utils::JSON::WriteHDF5 ( const std::string& outputFile, const DataSet& dataset, int compressionLevel )
{
	return HDF5DataSetIO::write(outputFile, dataset, compressionLevel);
}

bool Utils::JSON::ReadHDF5 ( const std::string& inputFile, DataSet& dataset )
{
	return HDF5DataSetIO::read(inputFile, dataset);
}