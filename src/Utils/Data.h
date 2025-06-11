#pragma once

#include <string>
#include <vector>
#include <queue>
#include <nlohmann/json.hpp>
#include <ofGraphics.h>
#include <ofSoundBuffer.h>

#define DATA_CHANGE_CHECK_1

#define DATA_NUM_STATS 7

namespace Acorex {
namespace Utils {

enum class Axis : int {
	X = 0,
	Y = 1,
	Z = 2,
	COLOR = 3,
	NONE = 4,
	MULTIPLE = 5
};

struct AudioData {
	std::vector<bool> loaded; // [file]
	std::vector<ofSoundBuffer> raw; // [file]
};

struct TimeData {
	std::vector<std::vector<std::vector<double>>> raw; // [file][timepoint][dimension] (first dimension is always time)
};

struct StatsData {
	std::vector<std::vector<std::vector<double>>> raw; // [file][dimension][statistic] (mean, stdDev, skewness, kurtosis, loPercent, midPercent, hiPercent)
	std::vector<std::vector<double>> reduced; // [file][dimension]
};

struct AnalysisSettings {
	int currentDimensionCount = 0;
	bool hasBeenReduced = false;
	bool bTime = false;
	bool bPitch = false;
	bool bLoudness = false;
	bool bShape = false;
	bool bMFCC = false;
	bool bTransport = false;
	int windowFFTSize = 1024;
	int hopFraction = 2;
	int nBands = 40;
	int nCoefs = 13;
	int minFreq = 20;
	int maxFreq = 5000;
	int sampleRate = 44100;
};

struct ReductionSettings {
	int dimensionReductionTarget = 3;
	int maxIterations = 200;
};

struct DataSet {
	int currentPointCount = 0;

	std::vector<std::string> dimensionNames; // [dimension]
	std::vector<std::string> statisticNames = { "Mean Average", "Standard Deviation", "Skewness", "Kurtosis", "Low Quartile", "Median", "High Quartile" }; // [statistic]
	std::vector<std::string> fileList; // [file]

	AudioData audio;

	TimeData time;

	StatsData stats;
	
	TransportData transport;

	AnalysisSettings analysisSettings;
};

struct PointFT {
	size_t file = 0;
	size_t time = 0;
};

struct AudioPlayhead {
	AudioPlayhead ( size_t ID, size_t file, size_t sample ) : playheadID ( ID ), fileIndex ( file ), sampleIndex ( sample ) { }

	size_t playheadID = 0;

	size_t fileIndex = 0;
	size_t sampleIndex = 0;

	bool crossfading = false;
	size_t jumpFileIndex = 0;
	size_t jumpSampleIndex = 0;
	size_t crossfadeCurrentSample = 0;
	size_t crossfadeSampleLength = 0;
	
	std::queue<size_t> triggerSamplePoints;
};

struct VisualPlayhead {
	VisualPlayhead ( size_t ID, size_t file, size_t sample ) : playheadID ( ID ), fileIndex ( file ), sampleIndex ( sample ) { }

	bool highlight = false;

	size_t playheadID = 0;

	size_t fileIndex = 0;
	size_t sampleIndex = 0;

	float position[3] = { 0.0, 0.0, 0.0 };

	ofColor color = ofColor ( 255, 255, 255, 255 );
	ofRectangle panelRect = ofRectangle ( 0, 0, 0, 0 );
};

struct TransportFrame {
	std::vector<double> magnitude;
	std::vector<double> phase;
	std::vector<double> dH;
	
	TransportFrame() = default;
	
	TransportFrame(size_t frameSize) 
		: magnitude(frameSize, 0.0)
		, phase(frameSize, 0.0)
		, dH(frameSize, 0.0) {}
		
	TransportFrame(const std::vector<double>& mag, 
	               const std::vector<double>& ph, 
	               const std::vector<double>& deriv)
		: magnitude(mag)
		, phase(ph)
		, dH(deriv) {}
		
	TransportFrame(std::vector<double>&& mag, 
	               std::vector<double>&& ph, 
	               std::vector<double>&& deriv)
		: magnitude(std::move(mag))
		, phase(std::move(ph))
		, dH(std::move(deriv)) {}
	
	size_t frameSize() const {
		return magnitude.size();
	}
	
	bool isValid() const {
		return !magnitude.empty() && 
		       magnitude.size() == phase.size() && 
		       magnitude.size() == dH.size();
	}
	
	void clear() {
		magnitude.clear();
		phase.clear();
		dH.clear();
	}
	
	void resize(size_t newSize) {
		magnitude.resize(newSize, 0.0);
		phase.resize(newSize, 0.0);
		dH.resize(newSize, 0.0);
	}
};

struct TransportData {
	std::vector<std::vector<TransportFrame>> frames;
	
	TransportData() = default;
	
	void clear() {
		frames.clear();
	}
	
	size_t fileCount() const {
		return frames.size();
	}
	
	size_t frameCount(size_t fileIndex) const {
		if (fileIndex < frames.size()) {
			return frames[fileIndex].size();
		}
		return 0;
	}
	
	bool hasFile(size_t fileIndex) const {
		return fileIndex < frames.size() && !frames[fileIndex].empty();
	}
	
	TransportFrame* getFrame(size_t fileIndex, size_t frameIndex) {
		if (fileIndex < frames.size() && frameIndex < frames[fileIndex].size()) {
			return &frames[fileIndex][frameIndex];
		}
		return nullptr;
	}
	
	const TransportFrame* getFrame(size_t fileIndex, size_t frameIndex) const {
		if (fileIndex < frames.size() && frameIndex < frames[fileIndex].size()) {
			return &frames[fileIndex][frameIndex];
		}
		return nullptr;
	}
	
	void addFile() {
		frames.emplace_back();
	}
	
	void reserveFiles(size_t count) {
		frames.reserve(count);
	}
	
	void reserveFrames(size_t fileIndex, size_t count) {
		if (fileIndex < frames.size()) {
			frames[fileIndex].reserve(count);
		}
	}
};

} // namespace Utils
} // namespace Acorex