#pragma once

#include <ofSoundBuffer.h>
#include <vector>
#include <mutex>
#include <algorithms/public/KDTree.hpp>
#include <data/TensorTypes.hpp>
#include <memory>
#include "AudioTransportN.hpp"

namespace Acorex {
namespace Explorer {

class AudioMorphEngine {
public:
    AudioMorphEngine();
    ~AudioMorphEngine();

    // Initialise the engine with audio processing parameters
    void Initialise(int sampleRate = 44100, size_t bufferSize = 512, float overlap = 0.75f);

    // Generate the next buffer of audio. Expects an already-allocated buffer with the desired
    // number of frames and channels (typically stereo).
    void Process(ofSoundBuffer &outBuffer);

    // KD-tree integration ----------------------------------------------
    using KDTree = fluid::algorithm::KDTree;
    using KNNResult = KDTree::KNNResult;
    using BarycentricWeights = AudioTransportN::BarycentricWeights;

    // Attach a pre-built KD-tree (shared ownership for thread-safety)
    void SetKDTree(const std::shared_ptr<KDTree>& kdTree);

    // Access the KD-tree (const) – returns nullptr if not set
    const std::shared_ptr<KDTree>& GetKDTree() const { return mKDTree; }

    /**
     * Query the KD-tree for the k nearest frames to the provided playhead vector.
     * @param query        Feature vector representing current playhead position
     * @param k            Number of nearest neighbours to return (default 3)
     * @param radius       Optional search radius (0 == unlimited)
     * @return             KDTree::KNNResult (distances + id pointers). Empty if tree invalid.
     */
    KNNResult QueryKNearest(const fluid::RealVector& query, size_t k = 3, double radius = 0.0) const;

    /**
     * Convert KD-tree distances into barycentric weights.
     * If gaussian==false, uses inverse distance weighting with exponent param (default 2).
     * If gaussian==true, uses Gaussian kernel with sigma=param (param must be >0, default 0.2).
     */
    BarycentricWeights CalculateWeights(const KNNResult& knn, bool gaussian=false, double param=2.0) const;

    int GetSampleRate() const { return mSampleRate; }
    size_t GetBufferSize() const { return mBufferSize; }
    float GetOverlap() const { return mOverlap; }

    // Reset internal state (clears buffers and counters)
    void Reset();

private:
    int mSampleRate;
    size_t mBufferSize;
    float mOverlap; // 0-1 range (e.g. 0.75 == 75 % overlap)

    // Buffers used for overlap-add synthesis and circular buffering
    std::vector<float> mOverlapBuffer;   // Scratch buffer for overlapped region
    std::vector<float> mCircularBuffer;  // Main circular buffer (mono)

    // Timing / indexing -----------------------------
    size_t mWritePos = 0;   // Next write index in circular buffer
    size_t mReadPos  = 0;   // Next read index for output (may lead write during overlap)
    uint64_t mSampleCounter = 0; // Running total samples processed (for diagnostics)

    // KD-tree instance (shared with other components)
    std::shared_ptr<KDTree> mKDTree;

    // Thread-safety for real-time operation (mutable for const query methods)
    mutable std::mutex mProcessMutex;

    static void NormalizeWeights(BarycentricWeights& w);
};

} // namespace Explorer
} // namespace Acorex 