/*
The MIT License (MIT)

Copyright (c) 2024 Elowyn Fearne

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <unordered_map>
#include <list>
#include <memory>
#include <mutex>
#include <chrono>
#include <atomic>
#include <vector>
#include <functional>

namespace Acorex {
namespace Explorer {

/**
 * Dynamic cache management system with LRU eviction and adaptive sizing
 * Optimized for audio buffer caching with real-time constraints
 */
class CacheManager {
public:
    struct CacheStats {
        std::atomic<size_t> hits{0};
        std::atomic<size_t> misses{0};
        std::atomic<size_t> evictions{0};
        std::atomic<size_t> currentSize{0};
        std::atomic<size_t> maxSize{0};
        std::atomic<double> hitRate{0.0};
        std::chrono::steady_clock::time_point lastUpdate;
        
        void updateHitRate() {
            size_t total = hits + misses;
            if (total > 0) {
                hitRate = static_cast<double>(hits) / total;
            }
        }
    };
    
    struct AudioBufferKey {
        size_t fileIndex;
        size_t frameIndex;
        
        bool operator==(const AudioBufferKey& other) const {
            return fileIndex == other.fileIndex && frameIndex == other.frameIndex;
        }
    };
    
    struct AudioBufferKeyHash {
        size_t operator()(const AudioBufferKey& key) const {
            return std::hash<size_t>()(key.fileIndex) ^ (std::hash<size_t>()(key.frameIndex) << 1);
        }
    };
    
    using AudioBuffer = std::vector<float>;
    using AudioBufferPtr = std::shared_ptr<AudioBuffer>;
    
private:
    // LRU cache node
    struct CacheNode {
        AudioBufferKey key;
        AudioBufferPtr buffer;
        size_t size;
        std::chrono::steady_clock::time_point lastAccess;
        size_t accessCount{0};
    };
    
    using CacheList = std::list<CacheNode>;
    using CacheMap = std::unordered_map<AudioBufferKey, CacheList::iterator, AudioBufferKeyHash>;
    
    CacheList mCacheList;
    CacheMap mCacheMap;
    mutable std::mutex mCacheMutex;
    
    CacheStats mStats;
    
    // Dynamic cache sizing parameters
    size_t mMinCacheSize = 64 * 1024 * 1024;     // 64MB minimum
    size_t mMaxCacheSize = 512 * 1024 * 1024;    // 512MB maximum
    size_t mCurrentMaxSize;
    
    // Performance monitoring
    std::chrono::steady_clock::time_point mLastResizeCheck;
    const std::chrono::seconds mResizeInterval{10};
    
    // Memory pressure callback
    std::function<size_t()> mAvailableMemoryCallback;
    
public:
    CacheManager(size_t initialMaxSize = 256 * 1024 * 1024);
    ~CacheManager();
    
    /**
     * Get audio buffer from cache or load it
     * @param key The buffer key (file + frame index)
     * @param loader Callback to load buffer if not in cache
     * @return Shared pointer to audio buffer
     */
    AudioBufferPtr get(const AudioBufferKey& key, 
                      std::function<AudioBufferPtr()> loader);
    
    /**
     * Prefetch audio buffers in background
     * @param keys List of buffer keys to prefetch
     * @param loader Callback to load each buffer
     */
    void prefetch(const std::vector<AudioBufferKey>& keys,
                 std::function<AudioBufferPtr(const AudioBufferKey&)> loader);
    
    /**
     * Clear the entire cache
     */
    void clear();
    
    /**
     * Get cache statistics
     */
    const CacheStats& getStats() const { return mStats; }
    
    /**
     * Set callback to query available system memory
     */
    void setAvailableMemoryCallback(std::function<size_t()> callback) {
        mAvailableMemoryCallback = callback;
    }
    
    /**
     * Manually trigger cache size adjustment
     */
    void adjustCacheSize();
    
    /**
     * Set cache size limits
     */
    void setSizeLimits(size_t minSize, size_t maxSize);
    
private:
    /**
     * Evict least recently used items to make space
     * @param requiredSpace Space needed in bytes
     */
    void evictLRU(size_t requiredSpace);
    
    /**
     * Move item to front of LRU list
     */
    void moveToFront(CacheList::iterator it);
    
    /**
     * Check and adjust cache size based on system memory and hit rate
     */
    void checkAndAdjustSize();
    
    /**
     * Calculate optimal cache size based on metrics
     */
    size_t calculateOptimalSize();
};

/**
 * Global cache manager singleton for audio buffers
 */
class GlobalCacheManager {
public:
    static CacheManager& getInstance() {
        static CacheManager instance;
        return instance;
    }
    
private:
    GlobalCacheManager() = default;
    ~GlobalCacheManager() = default;
    GlobalCacheManager(const GlobalCacheManager&) = delete;
    GlobalCacheManager& operator=(const GlobalCacheManager&) = delete;
};

} // namespace Explorer
} // namespace Acorex