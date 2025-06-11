#pragma once

#include <unordered_map>
#include <list>
#include <memory>
#include <shared_mutex>
#include <atomic>
#include <chrono>
#include <functional>
#include <optional>
#include "Data.h"

namespace Acorex {
namespace Utils {

// Cache statistics structure
struct CacheStats {
    std::atomic<size_t> hits{0};
    std::atomic<size_t> misses{0};
    std::atomic<size_t> evictions{0};
    std::atomic<size_t> currentSize{0};
    size_t maxSize{0};
    
    double getHitRate() const {
        size_t total = hits.load() + misses.load();
        return total > 0 ? static_cast<double>(hits.load()) / total : 0.0;
    }
    
    void reset() {
        hits.store(0);
        misses.store(0);
        evictions.store(0);
    }
};

// Generic LRU cache implementation
template<typename Key, typename Value>
class LRUCache {
public:
    using KeyType = Key;
    using ValueType = Value;
    using Clock = std::chrono::steady_clock;
    using TimePoint = Clock::time_point;
    
    struct CacheEntry {
        Key key;
        Value value;
        TimePoint lastAccess;
        
        CacheEntry(const Key& k, const Value& v) 
            : key(k), value(v), lastAccess(Clock::now()) {}
    };
    
    using CacheList = std::list<CacheEntry>;
    using CacheMap = std::unordered_map<Key, typename CacheList::iterator>;
    
    explicit LRUCache(size_t maxSize) : mMaxSize(maxSize) {
        if (maxSize == 0) {
            throw std::invalid_argument("Cache size must be greater than 0");
        }
    }
    
    // Get value from cache, returns empty optional if not found
    std::optional<Value> get(const Key& key) {
        auto it = mCacheMap.find(key);
        if (it == mCacheMap.end()) {
            return std::nullopt;
        }
        
        // Move to front (most recently used)
        mCacheList.splice(mCacheList.begin(), mCacheList, it->second);
        it->second->lastAccess = Clock::now();
        
        return it->second->value;
    }
    
    // Put value into cache
    void put(const Key& key, const Value& value) {
        auto it = mCacheMap.find(key);
        
        if (it != mCacheMap.end()) {
            // Update existing entry
            it->second->value = value;
            it->second->lastAccess = Clock::now();
            mCacheList.splice(mCacheList.begin(), mCacheList, it->second);
        } else {
            // Add new entry
            if (mCacheList.size() >= mMaxSize) {
                evictLRU();
            }
            
            mCacheList.emplace_front(key, value);
            mCacheMap[key] = mCacheList.begin();
        }
    }
    
    // Remove specific key from cache
    bool remove(const Key& key) {
        auto it = mCacheMap.find(key);
        if (it == mCacheMap.end()) {
            return false;
        }
        
        mCacheList.erase(it->second);
        mCacheMap.erase(it);
        return true;
    }
    
    // Clear all entries
    void clear() {
        mCacheList.clear();
        mCacheMap.clear();
    }
    
    // Get current cache size
    size_t size() const {
        return mCacheList.size();
    }
    
    // Check if cache contains key
    bool contains(const Key& key) const {
        return mCacheMap.find(key) != mCacheMap.end();
    }
    
    // Get maximum cache size
    size_t maxSize() const {
        return mMaxSize;
    }
    
    // Set new maximum size (may trigger evictions)
    void setMaxSize(size_t newSize) {
        if (newSize == 0) {
            throw std::invalid_argument("Cache size must be greater than 0");
        }
        
        mMaxSize = newSize;
        while (mCacheList.size() > mMaxSize) {
            evictLRU();
        }
    }
    
protected:
    // Evict least recently used entry
    void evictLRU() {
        if (!mCacheList.empty()) {
            auto last = mCacheList.end();
            --last;
            mCacheMap.erase(last->key);
            mCacheList.pop_back();
        }
    }
    
private:
    size_t mMaxSize;
    CacheList mCacheList;  // Front = most recently used, Back = least recently used
    CacheMap mCacheMap;
};

// Thread-safe LRU cache wrapper
template<typename Key, typename Value>
class ThreadSafeLRUCache {
public:
    using KeyType = Key;
    using ValueType = Value;
    using MissHandler = std::function<std::optional<Value>(const Key&)>;
    
    explicit ThreadSafeLRUCache(size_t maxSize) 
        : mCache(maxSize), mStats{} {
        mStats.maxSize = maxSize;
    }
    
    // Get value with read lock
    std::optional<Value> get(const Key& key) {
        {
            std::shared_lock<std::shared_mutex> lock(mMutex);
            auto result = mCache.get(key);
            if (result.has_value()) {
                mStats.hits.fetch_add(1);
                return result;
            }
        }
        
        // Cache miss
        mStats.misses.fetch_add(1);
        
        // Try miss handler if set
        if (mMissHandler) {
            auto value = mMissHandler(key);
            if (value.has_value()) {
                put(key, value.value());
                return value;
            }
        }
        
        return std::nullopt;
    }
    
    // Put value with write lock
    void put(const Key& key, const Value& value) {
        size_t oldSize, newSize;
        {
            std::unique_lock<std::shared_mutex> lock(mMutex);
            oldSize = mCache.size();
            mCache.put(key, value);
            newSize = mCache.size();
        }
        
        // Update stats
        if (newSize < oldSize) {
            mStats.evictions.fetch_add(oldSize - newSize);
        }
        mStats.currentSize.store(newSize);
    }
    
    // Remove with write lock
    bool remove(const Key& key) {
        bool removed;
        size_t newSize;
        {
            std::unique_lock<std::shared_mutex> lock(mMutex);
            removed = mCache.remove(key);
            newSize = mCache.size();
        }
        mStats.currentSize.store(newSize);
        return removed;
    }
    
    // Clear with write lock
    void clear() {
        {
            std::unique_lock<std::shared_mutex> lock(mMutex);
            mCache.clear();
        }
        mStats.currentSize.store(0);
    }
    
    // Thread-safe size query
    size_t size() const {
        std::shared_lock<std::shared_mutex> lock(mMutex);
        return mCache.size();
    }
    
    // Thread-safe contains check
    bool contains(const Key& key) const {
        std::shared_lock<std::shared_mutex> lock(mMutex);
        return mCache.contains(key);
    }
    
    // Set new max size
    void setMaxSize(size_t newSize) {
        size_t evicted = 0;
        {
            std::unique_lock<std::shared_mutex> lock(mMutex);
            size_t oldCacheSize = mCache.size();
            mCache.setMaxSize(newSize);
            size_t newCacheSize = mCache.size();
            if (newCacheSize < oldCacheSize) {
                evicted = oldCacheSize - newCacheSize;
            }
        }
        
        mStats.maxSize = newSize;
        if (evicted > 0) {
            mStats.evictions.fetch_add(evicted);
            mStats.currentSize.store(size());
        }
    }
    
    // Set miss handler
    void setMissHandler(MissHandler handler) {
        std::unique_lock<std::shared_mutex> lock(mMutex);
        mMissHandler = handler;
    }
    
    // Get statistics
    CacheStats getStats() const {
        return mStats;
    }
    
    // Reset statistics
    void resetStats() {
        mStats.reset();
    }
    
private:
    mutable std::shared_mutex mMutex;
    LRUCache<Key, Value> mCache;
    CacheStats mStats;
    MissHandler mMissHandler;
};

// Frame key for cache lookup
struct FrameKey {
    size_t fileIndex;
    size_t frameIndex;
    
    bool operator==(const FrameKey& other) const {
        return fileIndex == other.fileIndex && frameIndex == other.frameIndex;
    }
};

// Hash function for FrameKey
struct FrameKeyHash {
    std::size_t operator()(const FrameKey& key) const {
        // Combine hashes of both indices
        std::size_t h1 = std::hash<size_t>{}(key.fileIndex);
        std::size_t h2 = std::hash<size_t>{}(key.frameIndex);
        return h1 ^ (h2 << 1);
    }
};

// Forward declarations
class TransportDataLoader;
class FrameMemoryPool;

// TransportFrame-specific cache
class TransportFrameCache {
public:
    static constexpr size_t DEFAULT_CACHE_SIZE = 80000;  // ~30min corpus at typical settings
    static constexpr size_t DEFAULT_POOL_SIZE = 100000;  // Memory pool size
    
    using FramePtr = std::shared_ptr<TransportFrame>;
    using Cache = ThreadSafeLRUCache<FrameKey, FramePtr>;
    
    explicit TransportFrameCache(size_t maxFrames = DEFAULT_CACHE_SIZE);
    ~TransportFrameCache();
    
    // Get frame from cache (may trigger loading)
    FramePtr getFrame(size_t fileIndex, size_t frameIndex);
    
    // Batch get for efficiency
    std::vector<FramePtr> getFrames(const std::vector<FrameKey>& keys);
    
    // Preload frames
    void preloadFrames(const std::vector<FrameKey>& keys);
    
    // Cache management
    void clear();
    void setMaxSize(size_t maxFrames);
    size_t size() const { return mCache.size(); }
    
    // Set data loader for cache misses
    void setDataLoader(std::shared_ptr<TransportDataLoader> loader);
    
    // Statistics
    CacheStats getStats() const { return mCache.getStats(); }
    void resetStats() { mCache.resetStats(); }
    
    // Memory pool stats
    struct MemoryPoolStats {
        size_t totalAllocated;
        size_t currentlyUsed;
        size_t peakUsed;
        size_t poolSize;
    };
    MemoryPoolStats getMemoryPoolStats() const;
    
    // Performance monitoring
    struct PerformanceMetrics {
        double avgLoadTimeMs;
        double avgAccessTimeUs;
        size_t totalLoads;
        size_t totalAccesses;
    };
    PerformanceMetrics getPerformanceMetrics() const;
    
private:
    Cache mCache;
    std::shared_ptr<TransportDataLoader> mDataLoader;
    std::unique_ptr<FrameMemoryPool> mMemoryPool;
    
    // Performance tracking
    mutable std::atomic<uint64_t> mTotalLoadTimeUs{0};
    mutable std::atomic<size_t> mTotalLoads{0};
    mutable std::atomic<uint64_t> mTotalAccessTimeNs{0};
    mutable std::atomic<size_t> mTotalAccesses{0};
    
    // Frame loading implementation
    FramePtr loadFrame(const FrameKey& key);
};

// Interface for loading TransportFrame data
class TransportDataLoader {
public:
    virtual ~TransportDataLoader() = default;
    
    // Load a single frame
    virtual std::shared_ptr<TransportFrame> loadFrame(size_t fileIndex, size_t frameIndex) = 0;
    
    // Batch load for efficiency
    virtual std::vector<std::shared_ptr<TransportFrame>> loadFrames(
        const std::vector<FrameKey>& keys) {
        std::vector<std::shared_ptr<TransportFrame>> frames;
        frames.reserve(keys.size());
        for (const auto& key : keys) {
            frames.push_back(loadFrame(key.fileIndex, key.frameIndex));
        }
        return frames;
    }
    
    // Check if frame exists
    virtual bool hasFrame(size_t fileIndex, size_t frameIndex) const = 0;
};

} // namespace Utils
} // namespace Acorex