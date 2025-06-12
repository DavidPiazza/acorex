#include "./AudioMorphEngine.h"
#include "ofLog.h"
#include <cmath>
#include <utility> // for std::move if needed
#include <numeric>
#include <limits>
#include <algorithm>
#include <data/FluidMemory.hpp>

using namespace Acorex;

namespace Acorex {
namespace Explorer {

AudioMorphEngine::AudioMorphEngine()
: mSampleRate(44100), mBufferSize(512), mOverlap(0.75f) {
    // Create default allocator
    mAllocator = new fluid::Allocator();
    
    // Initialize with default FFT parameters
    mWindowSize = 2048;
    mFFTSize = 2048;
    mHopSize = static_cast<size_t>(mWindowSize * (1.0f - mOverlap));  // 512 for 75% overlap
    
    // Create AudioTransportN instance
    mTransportN = std::make_unique<AudioTransportN>(mFFTSize, *mAllocator);
}

AudioMorphEngine::~AudioMorphEngine() {
    if (mAllocator) {
        delete mAllocator;
        mAllocator = nullptr;
    }
}

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
    
    // Calculate hop size from overlap percentage
    mHopSize = static_cast<size_t>(mWindowSize * (1.0f - overlap));
    
    // Initialize AudioTransportN with our FFT parameters
    mTransportN->initN(mWindowSize, mFFTSize, mHopSize);
    
    // Size circular buffer to minimize latency while avoiding underruns
    // For <10ms latency at 44.1kHz, we need < 441 samples buffered
    // Use bufferSize + small margin for jitter handling
    size_t minBufferSize = bufferSize + 64; // Small margin for timing jitter
    size_t targetBufferSize = static_cast<size_t>(mSampleRate * 0.008); // 8ms target
    size_t circularBufferSize = std::max(minBufferSize, targetBufferSize);
    mOutputBuffer.resize(circularBufferSize);
    
    // Initialize overlap accumulator
    mOverlapAccumulator.assign(mWindowSize + mHopSize, 0.0f);
    
    // Initialize window function
    InitializeWindows();
    
    // Reset frame state
    mFrameState.currentFrame = 0;
    mFrameState.samplesUntilNextFrame = 0;
    mFrameState.lastProcessTime = std::chrono::steady_clock::now();
    mFrameState.averageLatency = 0.0;
    mFrameState.latencyMeasureCount = 0;
    
    // Reset performance counters
    mCurrentLatencyMs = 0.0;
    mUnderrunCount = 0;
    mOverrunCount = 0;

    ofLogNotice("AudioMorphEngine") << "Initialised with sampleRate=" << sampleRate
                                    << " bufferSize=" << bufferSize
                                    << " overlap=" << overlap * 100.0f << "%"
                                    << " windowSize=" << mWindowSize
                                    << " hopSize=" << mHopSize;
}

void AudioMorphEngine::Reset() {
    std::lock_guard<std::mutex> lock(mProcessMutex);
    
    // Reset circular buffer
    mOutputBuffer.reset();
    
    // Clear overlap accumulator
    std::fill(mOverlapAccumulator.begin(), mOverlapAccumulator.end(), 0.0f);
    
    // Reset frame state
    mFrameState.currentFrame = 0;
    mFrameState.samplesUntilNextFrame = 0;
    mFrameState.lastProcessTime = std::chrono::steady_clock::now();
    
    // Reset performance metrics
    mCurrentLatencyMs = 0.0;
    mUnderrunCount = 0;
    mOverrunCount = 0;
}

void AudioMorphEngine::Process(ofSoundBuffer &outBuffer) {
    std::lock_guard<std::mutex> lock(mProcessMutex);

    // Ensure the buffer is sized correctly
    if (static_cast<size_t>(outBuffer.getNumFrames()) != mBufferSize) {
        ofLogWarning("AudioMorphEngine") << "Provided buffer size (" << outBuffer.getNumFrames()
                                           << ") does not match engine buffer size (" << mBufferSize
                                           << "). Resizing (clears data).";
        outBuffer.setSampleRate(mSampleRate);
        outBuffer.resize(mBufferSize, outBuffer.getNumChannels());
    }
    
    size_t framesToProcess = outBuffer.getNumFrames();
    
    // Calculate minimum buffer threshold to avoid underruns
    size_t minBufferThreshold = framesToProcess + mHopSize;
    
    // Process frames proactively to maintain buffer
    while (mOutputBuffer.availableSamples() < minBufferThreshold) {
        // Check if it's time to process next frame
        if (mFrameState.samplesUntilNextFrame == 0) {
            ProcessNextFrame();
            mFrameState.samplesUntilNextFrame = mHopSize;
        }
        
        // If we've already processed a frame, wait for next hop
        if (mFrameState.samplesUntilNextFrame > 0) {
            break;
        }
    }
    
    // Read from circular buffer to output
    std::vector<float> monoOutput(framesToProcess);
    size_t samplesRead = 0;
    
    if (mOutputBuffer.availableSamples() >= framesToProcess) {
        mOutputBuffer.read(monoOutput.data(), framesToProcess);
        samplesRead = framesToProcess;
    } else {
        // Underrun - read what we can and pad with zeros
        size_t available = mOutputBuffer.availableSamples();
        if (available > 0) {
            mOutputBuffer.read(monoOutput.data(), available);
        }
        std::fill(monoOutput.begin() + available, monoOutput.end(), 0.0f);
        mUnderrunCount++;
        ofLogWarning("AudioMorphEngine") << "Buffer underrun: only " << available 
                                         << " samples available, needed " << framesToProcess;
    }
    
    // Copy mono output to all channels
    for (size_t frame = 0; frame < framesToProcess; ++frame) {
        for (size_t ch = 0; ch < outBuffer.getNumChannels(); ++ch) {
            outBuffer.getSample(frame, ch) = monoOutput[frame];
        }
    }
    
    // Update latency metrics
    UpdateLatencyMetrics();
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

// Circular Buffer Implementation
void AudioMorphEngine::CircularBuffer::resize(size_t newCapacity) {
    data.resize(newCapacity);
    capacity = newCapacity;
    reset();
}

size_t AudioMorphEngine::CircularBuffer::availableSamples() const {
    size_t write = writePos.load();
    size_t read = readPos.load();
    
    if (write >= read) {
        return write - read;
    } else {
        return capacity - read + write;
    }
}

size_t AudioMorphEngine::CircularBuffer::freeSpace() const {
    return capacity - availableSamples() - 1; // -1 to distinguish full from empty
}

bool AudioMorphEngine::CircularBuffer::write(const float* samples, size_t count) {
    if (count > freeSpace()) {
        return false; // Would overrun
    }
    
    size_t write = writePos.load();
    
    for (size_t i = 0; i < count; ++i) {
        data[write] = samples[i];
        write = (write + 1) % capacity;
    }
    
    writePos.store(write);
    return true;
}

bool AudioMorphEngine::CircularBuffer::read(float* samples, size_t count) {
    if (count > availableSamples()) {
        return false; // Would underrun
    }
    
    size_t read = readPos.load();
    
    for (size_t i = 0; i < count; ++i) {
        samples[i] = data[read];
        read = (read + 1) % capacity;
    }
    
    readPos.store(read);
    return true;
}

void AudioMorphEngine::CircularBuffer::reset() {
    writePos.store(0);
    readPos.store(0);
    std::fill(data.begin(), data.end(), 0.0f);
}

// Window function initialization
void AudioMorphEngine::InitializeWindows() {
    mWindowFunction.resize(mWindowSize);
    
    // Generate Hann window for overlap-add synthesis
    for (size_t i = 0; i < mWindowSize; ++i) {
        double phase = 2.0 * M_PI * i / (mWindowSize - 1);
        mWindowFunction[i] = 0.5f * (1.0f - std::cos(phase));
    }
    
    // For STFT with overlap-add, we typically don't need to normalize the window
    // The analysis window (in STFT) and synthesis window (in ISTFT) handle this
    // Our window is used for smooth amplitude enveloping only
    
    // Optional: Apply sqrt for amplitude preservation in overlap-add
    // This ensures that overlapped windows sum to approximately 1.0
    for (float& w : mWindowFunction) {
        w = std::sqrt(w);
    }
}

// Apply window function to input frame
void AudioMorphEngine::ApplyWindow(const float* input, float* output, size_t size) {
    size_t windowSize = std::min(size, mWindowFunction.size());
    for (size_t i = 0; i < windowSize; ++i) {
        output[i] = input[i] * mWindowFunction[i];
    }
    // Zero-pad if needed
    for (size_t i = windowSize; i < size; ++i) {
        output[i] = 0.0f;
    }
}

// Process next frame using AudioTransportN
void AudioMorphEngine::ProcessNextFrame() {
    // Record frame processing start time
    auto frameStartTime = std::chrono::steady_clock::now();
    
    // TODO: Get input frames from KD-tree query results
    // For now, we'll generate a test signal
    
    if (!mKDTree || !mTransportN->initialized()) {
        // No KD-tree or transport not initialized - generate silence
        std::vector<float> silentFrame(mWindowSize, 0.0f);
        AccumulateOverlapAdd(silentFrame.data(), mWindowSize);
        return;
    }
    
    // Query KD-tree for nearest neighbors
    // TODO: Get query point from playhead position
    fluid::RealVector queryPoint(mKDTree->dims());
    // Initialize to zeros
    for (int i = 0; i < queryPoint.size(); ++i) {
        queryPoint[i] = 0.0;
    }
    
    auto knnResult = QueryKNearest(queryPoint, mNumNeighbors); // Get configured number of nearest neighbors
    if (knnResult.first.size() == 0) {
        // No results - generate silence
        std::vector<float> silentFrame(mWindowSize, 0.0f);
        AccumulateOverlapAdd(silentFrame.data(), mWindowSize);
        return;
    }
    
    // Calculate barycentric weights from distances
    auto weights = CalculateWeights(knnResult, false, 2.0);
    
    // TODO: Load actual audio frames from corpus using knnResult indices
    // For now, generate test frames
    std::vector<fluid::RealVectorView> inputFrames;
    std::vector<std::vector<double>> frameData(knnResult.first.size());
    
    for (size_t i = 0; i < knnResult.first.size(); ++i) {
        frameData[i].resize(mWindowSize, 0.0);
        // TODO: Load real audio data
        // For testing, generate different frequencies
        double freq = 440.0 * (i + 1);
        for (size_t j = 0; j < mWindowSize; ++j) {
            frameData[i][j] = 0.1 * std::sin(2.0 * M_PI * freq * j / mSampleRate);
        }
        inputFrames.emplace_back(frameData[i].data(), 0, mWindowSize);
    }
    
    // Process through AudioTransportN
    fluid::RealMatrix outputMatrix(2, mWindowSize);
    mTransportN->processFrameN(inputFrames, weights, outputMatrix);
    
    // Extract audio output (first row) and window (second row)
    std::vector<float> morphedFrame(mWindowSize);
    for (fluid::index i = 0; i < mWindowSize; ++i) {
        morphedFrame[i] = outputMatrix(0, i);
    }
    
    // Apply our window function and accumulate
    std::vector<float> windowedFrame(mWindowSize);
    ApplyWindow(morphedFrame.data(), windowedFrame.data(), mWindowSize);
    AccumulateOverlapAdd(windowedFrame.data(), mWindowSize);
    
    // Update frame counter
    mFrameState.currentFrame++;
    
    // Calculate frame processing latency
    auto frameEndTime = std::chrono::steady_clock::now();
    auto frameDuration = std::chrono::duration_cast<std::chrono::microseconds>
                        (frameEndTime - frameStartTime).count() / 1000.0;
    
    // Update running average latency
    mFrameState.averageLatency = (mFrameState.averageLatency * mFrameState.latencyMeasureCount + frameDuration) 
                                / (mFrameState.latencyMeasureCount + 1);
    mFrameState.latencyMeasureCount++;
}

// Accumulate windowed frame into overlap buffer
void AudioMorphEngine::AccumulateOverlapAdd(const float* frame, size_t frameSize) {
    // Shift overlap accumulator by hop size
    std::copy(mOverlapAccumulator.begin() + mHopSize, 
              mOverlapAccumulator.end(), 
              mOverlapAccumulator.begin());
    
    // Clear the tail
    std::fill(mOverlapAccumulator.end() - mHopSize, 
              mOverlapAccumulator.end(), 
              0.0f);
    
    // Add new frame with overlap
    for (size_t i = 0; i < frameSize && i < mOverlapAccumulator.size(); ++i) {
        mOverlapAccumulator[i] += frame[i];
    }
    
    // Extract hop size samples to circular buffer
    if (!mOutputBuffer.write(mOverlapAccumulator.data(), mHopSize)) {
        ofLogWarning("AudioMorphEngine") << "Failed to write overlap-add output to circular buffer";
    }
}

// Update latency metrics
void AudioMorphEngine::UpdateLatencyMetrics() {
    // Calculate current buffer latency
    size_t bufferedSamples = mOutputBuffer.availableSamples();
    double bufferLatencyMs = (bufferedSamples * 1000.0) / mSampleRate;
    
    // Add processing latency
    double totalLatencyMs = bufferLatencyMs + mFrameState.averageLatency;
    
    // Update atomic latency value
    mCurrentLatencyMs.store(totalLatencyMs);
    
    // Log if latency exceeds target
    if (totalLatencyMs > 10.0 && mFrameState.latencyMeasureCount % 100 == 0) {
        ofLogWarning("AudioMorphEngine") << "Latency exceeds target: " 
                                         << totalLatencyMs << "ms (target: <10ms)";
    }
}

// Set FFT parameters
void AudioMorphEngine::SetFFTParameters(size_t windowSize, size_t fftSize, size_t hopSize) {
    std::lock_guard<std::mutex> lock(mProcessMutex);
    
    // Validate parameters
    if (windowSize == 0 || fftSize < windowSize || hopSize == 0 || hopSize > windowSize) {
        ofLogError("AudioMorphEngine") << "Invalid FFT parameters";
        return;
    }
    
    mWindowSize = windowSize;
    mFFTSize = fftSize;
    mHopSize = hopSize;
    mOverlap = 1.0f - (float)hopSize / windowSize;
    
    // Reinitialize if already initialized
    if (mTransportN && mTransportN->initialized()) {
        mTransportN->initN(mWindowSize, mFFTSize, mHopSize);
        InitializeWindows();
        mOverlapAccumulator.resize(mWindowSize + mHopSize);
        std::fill(mOverlapAccumulator.begin(), mOverlapAccumulator.end(), 0.0f);
    }
    
    ofLogNotice("AudioMorphEngine") << "FFT parameters updated: windowSize=" << windowSize
                                    << " fftSize=" << fftSize 
                                    << " hopSize=" << hopSize
                                    << " (" << mOverlap * 100 << "% overlap)";
}

// -------------------------------------------------------------------

} // namespace Explorer
} // namespace Acorex 