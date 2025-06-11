#pragma once

#include <ofSoundBuffer.h>
#include <vector>
#include <mutex>

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

    // Placeholder setters/getters for the KD-tree that will be integrated in later subtasks
    void *GetKDTree() const { return mKDTree; }
    void SetKDTree(void *kdTree) { mKDTree = kdTree; }

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

    // KD-tree pointer – real type will be introduced in a later subtask
    void *mKDTree = nullptr;

    // Thread-safety for real-time operation
    std::mutex mProcessMutex;
};

} // namespace Explorer
} // namespace Acorex 