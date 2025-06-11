#include "./AudioMorphEngine.h"
#include "ofLog.h"
#include <cmath>

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

} // namespace Explorer
} // namespace Acorex 