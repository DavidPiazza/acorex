#include "./AudioMorphEngine.h"
#include "ofLog.h"
#include <cmath>
#include <utility> // for std::move if needed
#include <numeric>
#include <limits>

using namespace Acorex;

namespace Acorex {
namespace Explorer {

AudioMorphEngine::AudioMorphEngine()
: mSampleRate(44100), mBufferSize(512), mOverlap(0.75f) {
    // Nothing else to do here yet
}

AudioMorphEngine::~AudioMorphEngine() = default;

void AudioMorphEngine::Initialise(int sampleRate, size_t bufferSize, float overlap) {
    std::lock_guard<std::mutex> lock(mProcessMutex);

    // --- Parameter validation ----------------------
    if (sampleRate <= 0) {
        ofLogWarning("AudioMorphEngine") << "Invalid sampleRate=" << sampleRate << ", falling back to 44100.";
        sampleRate = 44100;
    }
    if (bufferSize == 0) {
        ofLogWarning("AudioMorphEngine") << "bufferSize may not be zero – defaulting to 512.";
        bufferSize = 512;
    }
    if (overlap <= 0.0f || overlap >= 1.0f) {
        ofLogWarning("AudioMorphEngine") << "overlap must be between 0 and 1 (exclusive). Clamping to 0.75.";
        overlap = 0.75f;
    }

    mSampleRate = sampleRate;
    mBufferSize = bufferSize;
    mOverlap = overlap;

    // Allocate buffers sized relative to the requested bufferSize
    mOverlapBuffer.assign(static_cast<size_t>(bufferSize * overlap), 0.0f);
    mCircularBuffer.assign(bufferSize, 0.0f);

    // Reset indices & counters
    mWritePos = 0;
    mReadPos = 0;
    mSampleCounter = 0;

    ofLogNotice("AudioMorphEngine") << "Initialised with sampleRate=" << sampleRate
                                      << " bufferSize=" << bufferSize
                                      << " overlap=" << overlap * 100.0f << "%";
}

void AudioMorphEngine::Reset() {
    std::lock_guard<std::mutex> lock(mProcessMutex);
    std::fill(mCircularBuffer.begin(), mCircularBuffer.end(), 0.0f);
    std::fill(mOverlapBuffer.begin(), mOverlapBuffer.end(), 0.0f);
    mWritePos = mReadPos = 0;
    mSampleCounter = 0;
}

void AudioMorphEngine::Process(ofSoundBuffer &outBuffer) {
    std::lock_guard<std::mutex> lock(mProcessMutex);

    // Ensure the buffer is sized correctly; resize if required (this should rarely happen)
    if (static_cast<size_t>(outBuffer.getNumFrames()) != mBufferSize) {
        ofLogWarning("AudioMorphEngine") << "Provided buffer size (" << outBuffer.getNumFrames()
                                           << ") does not match engine buffer size (" << mBufferSize
                                           << "). Resizing (clears data).";
        outBuffer.setSampleRate(mSampleRate);
        outBuffer.resize(mBufferSize, outBuffer.getNumChannels());
    }

    // Placeholder: For now, emit silence but update counters so that unit tests can
    // check timing-related values without implementing full synthesis yet.
    for (size_t frame = 0; frame < outBuffer.getNumFrames(); ++frame) {
        // Advance write position (simulate writing silence into circular buffer)
        mCircularBuffer[mWritePos] = 0.0f;
        mWritePos = (mWritePos + 1) % mCircularBuffer.size();

        // Read from circular buffer into all channels
        float sample = mCircularBuffer[mReadPos];
        mReadPos = (mReadPos + 1) % mCircularBuffer.size();

        for (size_t ch = 0; ch < outBuffer.getNumChannels(); ++ch) {
            outBuffer.getSample(frame, ch) = sample;
        }

        ++mSampleCounter;
    }
}

// ================= KD-tree Integration =============================

void Acorex::Explorer::AudioMorphEngine::SetKDTree(const std::shared_ptr<KDTree>& kdTree) {
    std::lock_guard<std::mutex> lock(mProcessMutex);
    mKDTree = kdTree;
    if (!mKDTree) {
        ofLogWarning("AudioMorphEngine") << "KDTree pointer set to nullptr. Querying will be disabled.";
    }
}

Acorex::Explorer::AudioMorphEngine::KNNResult
Acorex::Explorer::AudioMorphEngine::QueryKNearest(const fluid::RealVector& query,
                                                  size_t k,
                                                  double radius) const {
    // Fast path: early check without locking if KDTree is absent
    if (!mKDTree) {
        ofLogWarning("AudioMorphEngine") << "QueryKNearest called but KDTree not set.";
        return {};
    }

    std::lock_guard<std::mutex> lock(mProcessMutex);

    if (!mKDTree || !mKDTree->initialized()) {
        ofLogWarning("AudioMorphEngine") << "KDTree not initialized. Returning empty result.";
        return {};
    }

    // Guard against mismatched dimensionality
    if (static_cast<size_t>(query.size()) != static_cast<size_t>(mKDTree->dims())) {
        ofLogWarning("AudioMorphEngine") << "Query vector dimension " << query.size()
                                          << " does not match KDTree dims " << mKDTree->dims()
                                          << ". Returning empty result.";
        return {};
    }

    // Perform search (radius == 0 => unlimited)
    try {
        auto result = mKDTree->kNearest(query, static_cast<fluid::index>(k), radius);
        return result;
    } catch (const std::exception& e) {
        ofLogError("AudioMorphEngine") << "Exception during KDTree query: " << e.what();
        return {};
    }
}

// -------------------------------------------------------------------

Acorex::Explorer::AudioMorphEngine::BarycentricWeights
Acorex::Explorer::AudioMorphEngine::CalculateWeights(const KNNResult& knn,
                                                      bool gaussian,
                                                      double param) const {
    BarycentricWeights weights(knn.first.size());
    size_t n = knn.first.size();
    if (n == 0) return weights; // defaults to 1 on first weight via ctor

    // Edge case: if any distance is 0, return one-hot weight
    for (size_t i = 0; i < n; ++i) {
        if (knn.first[i] == 0.0) {
            weights.weights.assign(n, 0.0);
            weights[i] = 1.0;
            return weights;
        }
    }

    std::vector<double> raw(n, 0.0);
    if (gaussian) {
        double sigma = (param > 0.0) ? param : 0.2;
        double twoSigma2 = 2.0 * sigma * sigma;
        for (size_t i = 0; i < n; ++i) {
            raw[i] = std::exp(- (knn.first[i] * knn.first[i]) / twoSigma2);
        }
    } else {
        double p = (param > 0.0) ? param : 2.0;
        for (size_t i = 0; i < n; ++i) {
            raw[i] = 1.0 / std::pow(knn.first[i], p);
        }
    }

    weights.weights = std::move(raw);
    NormalizeWeights(weights);
    return weights;
}

void Acorex::Explorer::AudioMorphEngine::NormalizeWeights(BarycentricWeights& w) {
    double sum = std::accumulate(w.weights.begin(), w.weights.end(), 0.0);
    if (sum <= std::numeric_limits<double>::epsilon()) {
        if (!w.weights.empty()) { w.weights.assign(w.weights.size(), 0.0); w.weights[0] = 1.0; }
        return;
    }
    for (double& v : w.weights) v /= sum;
}

// -------------------------------------------------------------------

} // namespace Explorer
} // namespace Acorex 