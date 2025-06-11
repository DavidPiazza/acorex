#include "TransportFrameCache.h"
#include <algorithm>
#include <cstring>
#include <numeric>
#include <functional>

namespace Acorex {
namespace Utils {

// Memory pool for TransportFrame allocation
class FrameMemoryPool {
public:
    struct Block {
        alignas(TransportFrame) char data[sizeof(TransportFrame)];
        bool inUse = false;
    };
    
    explicit FrameMemoryPool(size_t poolSize) 
        : mPoolSize(poolSize), mBlocks(poolSize) {
    }
    
    std::shared_ptr<TransportFrame> allocate(size_t frameSize) {
        std::lock_guard<std::mutex> lock(mMutex);
        
        // Find free block
        auto it = std::find_if(mBlocks.begin(), mBlocks.end(), 
            [](const Block& b) { return !b.inUse; });
        
        if (it != mBlocks.end()) {
            it->inUse = true;
            mCurrentlyUsed++;
            mPeakUsed = std::max(mPeakUsed, mCurrentlyUsed);
            
            // Construct TransportFrame in-place
            TransportFrame* frame = new (it->data) TransportFrame(frameSize);
            
            // Return shared_ptr with custom deleter
            size_t index = std::distance(mBlocks.begin(), it);
            return std::shared_ptr<TransportFrame>(frame, 
                [this, index](TransportFrame* f) {
                    this->deallocate(f, index);
                });
        }
        
        // Pool exhausted, fall back to regular allocation
        mTotalAllocated++;
        return std::make_shared<TransportFrame>(frameSize);
    }
    
    void deallocate(TransportFrame* frame, size_t blockIndex) {
        frame->~TransportFrame();
        
        std::lock_guard<std::mutex> lock(mMutex);
        if (blockIndex < mBlocks.size()) {
            mBlocks[blockIndex].inUse = false;
            mCurrentlyUsed--;
        }
    }
    
    size_t getTotalAllocated() const {
        std::lock_guard<std::mutex> lock(mMutex);
        return mTotalAllocated;
    }
    
    size_t getCurrentlyUsed() const {
        std::lock_guard<std::mutex> lock(mMutex);
        return mCurrentlyUsed;
    }
    
    size_t getPeakUsed() const {
        std::lock_guard<std::mutex> lock(mMutex);
        return mPeakUsed;
    }
    
    size_t getPoolSize() const {
        return mPoolSize;
    }
    
private:
    mutable std::mutex mMutex;
    size_t mPoolSize;
    std::vector<Block> mBlocks;
    size_t mTotalAllocated = 0;
    size_t mCurrentlyUsed = 0;
    size_t mPeakUsed = 0;
};

// TransportFrameCache implementation
TransportFrameCache::TransportFrameCache(size_t maxFrames)
    : mCache(maxFrames)
    , mMemoryPool(std::make_unique<FrameMemoryPool>(DEFAULT_POOL_SIZE)) {
    
    // Set up miss handler
    mCache.setMissHandler([this](const FrameKey& key) -> std::optional<FramePtr> {
        auto frame = loadFrame(key);
        return frame ? std::optional<FramePtr>(frame) : std::nullopt;
    });
}

TransportFrameCache::~TransportFrameCache() = default;

TransportFrameCache::FramePtr TransportFrameCache::getFrame(size_t fileIndex, size_t frameIndex) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    FrameKey key{fileIndex, frameIndex};
    auto result = mCache.get(key);
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count();
    
    mTotalAccessTimeNs.fetch_add(duration);
    mTotalAccesses.fetch_add(1);
    
    return result.value_or(nullptr);
}

std::vector<TransportFrameCache::FramePtr> TransportFrameCache::getFrames(
    const std::vector<FrameKey>& keys) {
    
    std::vector<FramePtr> frames;
    frames.reserve(keys.size());
    
    for (const auto& key : keys) {
        frames.push_back(getFrame(key.fileIndex, key.frameIndex));
    }
    
    return frames;
}

void TransportFrameCache::preloadFrames(const std::vector<FrameKey>& keys) {
    if (!mDataLoader) return;
    
    // Filter out already cached frames
    std::vector<FrameKey> toLoad;
    toLoad.reserve(keys.size());
    
    for (const auto& key : keys) {
        if (!mCache.contains(key)) {
            toLoad.push_back(key);
        }
    }
    
    if (toLoad.empty()) return;
    
    // Batch load from data loader
    auto frames = mDataLoader->loadFrames(toLoad);
    
    // Add to cache
    for (size_t i = 0; i < toLoad.size() && i < frames.size(); ++i) {
        if (frames[i]) {
            mCache.put(toLoad[i], frames[i]);
        }
    }
}

void TransportFrameCache::clear() {
    mCache.clear();
}

void TransportFrameCache::setMaxSize(size_t maxFrames) {
    mCache.setMaxSize(maxFrames);
}

void TransportFrameCache::setDataLoader(std::shared_ptr<TransportDataLoader> loader) {
    mDataLoader = loader;
}

TransportFrameCache::MemoryPoolStats TransportFrameCache::getMemoryPoolStats() const {
    return {
        mMemoryPool->getTotalAllocated(),
        mMemoryPool->getCurrentlyUsed(),
        mMemoryPool->getPeakUsed(),
        mMemoryPool->getPoolSize()
    };
}

TransportFrameCache::PerformanceMetrics TransportFrameCache::getPerformanceMetrics() const {
    PerformanceMetrics metrics{};
    
    size_t loads = mTotalLoads.load();
    if (loads > 0) {
        metrics.avgLoadTimeMs = static_cast<double>(mTotalLoadTimeUs.load()) / (loads * 1000.0);
    }
    
    size_t accesses = mTotalAccesses.load();
    if (accesses > 0) {
        metrics.avgAccessTimeUs = static_cast<double>(mTotalAccessTimeNs.load()) / (accesses * 1000.0);
    }
    
    metrics.totalLoads = loads;
    metrics.totalAccesses = accesses;
    
    return metrics;
}

TransportFrameCache::FramePtr TransportFrameCache::loadFrame(const FrameKey& key) {
    if (!mDataLoader) {
        return nullptr;
    }
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    auto frame = mDataLoader->loadFrame(key.fileIndex, key.frameIndex);
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
    
    mTotalLoadTimeUs.fetch_add(duration);
    mTotalLoads.fetch_add(1);
    
    return frame;
}

} // namespace Utils
} // namespace Acorex