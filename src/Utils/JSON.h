#pragma once

#include "Utils/Data.h"
#include <nlohmann/json.hpp>

namespace Acorex {
namespace Utils {

class JSON {

#define TO_J( x ) {#x, a.x}
#define TO_J_SETTINGS( x ) {#x, a.analysisSettings.x}

#define TO_A( x ) j.at ( #x ).get_to ( a.x )
#define TO_A_SETTINGS( x ) j.at ( #x ).get_to ( a.analysisSettings.x )

public:
	JSON ( ) { };
	~JSON ( ) { };

	bool Write ( const std::string& outputFile, const DataSet& dataset );

	bool Read ( const std::string& inputFile, DataSet& dataset );
	bool Read ( const std::string& inputFile, AnalysisSettings& settings );
	
	// Hybrid I/O methods using memory-mapped files for transport data
	bool WriteHybrid ( const std::string& outputFile, const DataSet& dataset );
	bool ReadHybrid ( const std::string& inputFile, DataSet& dataset );
	
	// HDF5-based I/O methods with compression support
	bool WriteHDF5 ( const std::string& outputFile, const DataSet& dataset, int compressionLevel = 6 );
	bool ReadHDF5 ( const std::string& inputFile, DataSet& dataset );
};

void to_json ( nlohmann::json& j, const DataSet& a );
void from_json ( const nlohmann::json& j, DataSet& a );

void to_json ( nlohmann::json& j, const AnalysisSettings& a );
void from_json ( const nlohmann::json& j, AnalysisSettings& a );

void to_json ( nlohmann::json& j, const TransportFrame& a );
void from_json ( const nlohmann::json& j, TransportFrame& a );

void to_json ( nlohmann::json& j, const TransportData& a );
void from_json ( const nlohmann::json& j, TransportData& a );

} // namespace Utils
} // namespace Acorex