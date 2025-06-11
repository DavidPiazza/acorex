#pragma once

#include "TransportFrameCache.h"
#include "TransportDataLoaders.h"
#include "../Explorer/AudioTransportN.hpp"
#include <memory>

namespace Acorex {
namespace Utils {

/**
 * Example integration of TransportFrameCache with AudioTransportN
 * Shows how the cache would be used in the audio morphing pipeline
 */
class CachedAudioMorphEngine {
public:
    explicit CachedAudioMorphEngine(const std::string& corpusPath, size_t cacheSize = 80000)
        : mCache(std::make_unique<TransportFrameCache>(cacheSize)) {
        
        // Create appropriate loader based on file type
        auto loader = createTransportDataLoader(corpusPath);
        mCache->setDataLoader(loader);
    }
    
    /**
     * Process morphing with cached frames
     * @param frameKeys Array of frame keys to morph between
     * @param weights Barycentric weights for morphing
     * @param output Output audio buffer
     */
    void processMorph(
        const std::vector<FrameKey>& frameKeys,
        const Explorer::AudioTransportN::BarycentricWeights& weights,
        fluid::RealMatrixView output) {
        
        // Sanity check
        if (frameKeys.size() != weights.size()) {
            throw std::invalid_argument("Frame keys and weights must have same size");
        }
        
        // Get frames from cache (may trigger loading)
        auto frames = mCache->getFrames(frameKeys);
        
        // Convert to format expected by AudioTransportN
        std::vector<fluid::RealVectorView> frameViews;
        frameViews.reserve(frames.size());
        
        for (const auto& frame : frames) {
            if (frame && frame->isValid()) {
                // Here we would convert TransportFrame to the format needed by AudioTransportN
                // For now, using magnitude as example
                frameViews.emplace_back(
                    fluid::RealVectorView(frame->magnitude.data(), frame->magnitude.size())
                );
            }
        }
        
        // Perform N-way morphing
        if (mAudioTransport && !frameViews.empty()) {
            mAudioTransport->processFrameN(frameViews, weights, output);
        }
    }
    
    /**
     * Preload frames for upcoming playback
     * Call this ahead of time to ensure smooth playback
     */
    void preloadRegion(size_t fileIndex, size_t startFrame, size_t endFrame) {
        std::vector<FrameKey> keys;
        keys.reserve(endFrame - startFrame);
        
        for (size_t i = startFrame; i < endFrame; ++i) {
            keys.push_back({fileIndex, i});
        }
        
        mCache->preloadFrames(keys);
    }
    
    /**
     * Get cache statistics for monitoring
     */
    CacheStats getCacheStats() const {
        return mCache->getStats();
    }
    
    /**
     * Get performance metrics
     */
    TransportFrameCache::PerformanceMetrics getPerformanceMetrics() const {
        return mCache->getPerformanceMetrics();
    }
    
    /**
     * Clear cache (useful when switching corpora)
     */
    void clearCache() {
        mCache->clear();
    }
    
    /**
     * Set cache size
     */
    void setCacheSize(size_t maxFrames) {
        mCache->setMaxSize(maxFrames);
    }
    
private:
    std::unique_ptr<TransportFrameCache> mCache;
    std::unique_ptr<Explorer::AudioTransportN> mAudioTransport;
};

/**
 * Background frame prefetcher
 * Predicts and loads frames based on playhead movement
 */
class FramePrefetcher {
public:
    FramePrefetcher(std::shared_ptr<TransportFrameCache> cache)
        : mCache(cache) {}
    
    /**
     * Update playhead position and prefetch nearby frames
     */
    void updatePlayhead(size_t fileIndex, size_t frameIndex, double velocity) {
        // Calculate prefetch window based on velocity
        size_t prefetchDistance = static_cast<size_t>(
            std::abs(velocity) * mPrefetchTimeSec * mSampleRate / mHopSize
        );
        prefetchDistance = std::max(prefetchDistance, mMinPrefetchFrames);
        prefetchDistance = std::min(prefetchDistance, mMaxPrefetchFrames);
        
        // Determine direction
        if (velocity > 0) {
            // Moving forward
            prefetchFrames(fileIndex, frameIndex, frameIndex + prefetchDistance);
        } else if (velocity < 0) {
            // Moving backward
            size_t start = frameIndex > prefetchDistance ? frameIndex - prefetchDistance : 0;
            prefetchFrames(fileIndex, start, frameIndex);
        } else {
            // Stationary - prefetch surrounding area
            size_t start = frameIndex > prefetchDistance/2 ? frameIndex - prefetchDistance/2 : 0;
            prefetchFrames(fileIndex, start, frameIndex + prefetchDistance/2);
        }
    }
    
    /**
     * Set prefetch parameters
     */
    void setPrefetchParams(double timeSec, size_t minFrames, size_t maxFrames) {
        mPrefetchTimeSec = timeSec;
        mMinPrefetchFrames = minFrames;
        mMaxPrefetchFrames = maxFrames;
    }
    
private:
    std::shared_ptr<TransportFrameCache> mCache;
    double mPrefetchTimeSec = 2.0;  // Prefetch 2 seconds ahead
    size_t mMinPrefetchFrames = 100;
    size_t mMaxPrefetchFrames = 1000;
    size_t mSampleRate = 44100;
    size_t mHopSize = 512;
    
    void prefetchFrames(size_t fileIndex, size_t startFrame, size_t endFrame) {
        std::vector<FrameKey> keys;
        keys.reserve(endFrame - startFrame);
        
        for (size_t i = startFrame; i < endFrame; ++i) {
            keys.push_back({fileIndex, i});
        }
        
        // Could do this asynchronously in production
        mCache->preloadFrames(keys);
    }
};

} // namespace Utils
} // namespace Acorex