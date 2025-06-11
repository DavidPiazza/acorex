#pragma once

#include <ofSoundBuffer.h>
#include <vector>
#include <mutex>
#include <algorithms/public/KDTree.hpp>
#include <data/TensorTypes.hpp>
#include <memory>
#include <atomic>
#include <chrono>
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
    
    // Performance monitoring
    double GetCurrentLatencyMs() const { return mCurrentLatencyMs.load(); }
    size_t GetUnderrunCount() const { return mUnderrunCount.load(); }
    size_t GetOverrunCount() const { return mOverrunCount.load(); }
    
    // Configure FFT parameters
    void SetFFTParameters(size_t windowSize, size_t fftSize, size_t hopSize);

private:
    int mSampleRate;
    size_t mBufferSize;
    float mOverlap; // 0-1 range (e.g. 0.75 == 75 % overlap)

    // FFT and window parameters
    size_t mWindowSize = 2048;      // STFT window size
    size_t mHopSize = 512;          // Hop size (25% of window for 75% overlap)
    size_t mFFTSize = 2048;         // FFT size (typically same as window)
    
    // Circular buffer for overlap-add synthesis
    struct CircularBuffer {
        std::vector<float> data;
        std::atomic<size_t> writePos{0};
        std::atomic<size_t> readPos{0};
        size_t capacity;
        
        CircularBuffer() : capacity(0) {}
        void resize(size_t newCapacity);
        size_t availableSamples() const;
        size_t freeSpace() const;
        bool write(const float* samples, size_t count);
        bool read(float* samples, size_t count);
        void reset();
    };
    
    CircularBuffer mOutputBuffer;      // Main output circular buffer
    std::vector<float> mOverlapAccumulator;  // Accumulates overlapped windows
    std::vector<float> mWindowFunction;      // Pre-computed window function
    
    // Frame processing state
    struct FrameState {
        size_t currentFrame = 0;        // Current frame index
        size_t samplesUntilNextFrame = 0; // Samples until next frame processing
        std::chrono::steady_clock::time_point lastProcessTime;
        double averageLatency = 0.0;    // Running average latency in ms
        size_t latencyMeasureCount = 0;
    };
    FrameState mFrameState;
    
    // AudioTransportN instance for spectral morphing
    std::unique_ptr<AudioTransportN> mTransportN;
    fluid::Allocator* mAllocator = nullptr;
    
    // KD-tree instance (shared with other components)
    std::shared_ptr<KDTree> mKDTree;

    // Thread-safety for real-time operation (mutable for const query methods)
    mutable std::mutex mProcessMutex;
    
    // Performance metrics
    std::atomic<double> mCurrentLatencyMs{0.0};
    std::atomic<size_t> mUnderrunCount{0};
    std::atomic<size_t> mOverrunCount{0};

    static void NormalizeWeights(BarycentricWeights& w);
    
    // Internal processing methods
    void ProcessNextFrame();
    void ApplyWindow(const float* input, float* output, size_t size);
    void AccumulateOverlapAdd(const float* frame, size_t frameSize);
    void InitializeWindows();
    void UpdateLatencyMetrics();
};

} // namespace Explorer
} // namespace Acorex 