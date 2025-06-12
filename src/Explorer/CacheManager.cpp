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

#include "CacheManager.h"
#include <algorithm>
#include <thread>
#include <future>

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#include <mach/mach.h>
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#elif defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#elif defined(__linux__)
#include <sys/sysinfo.h>
#endif

namespace Acorex {
namespace Explorer {

// Helper function to get available system memory
static size_t getAvailableSystemMemory() {
#ifdef __APPLE__
    vm_size_t page_size;
    mach_port_t mach_port = mach_host_self();
    vm_statistics64_data_t vm_stat;
    mach_msg_type_number_t host_size = sizeof(vm_stat) / sizeof(natural_t);
    
    if (host_statistics64(mach_port, HOST_VM_INFO64, (host_info64_t)&vm_stat, &host_size) != KERN_SUCCESS) {
        return 0;
    }
    
    if (host_page_size(mach_port, &page_size) != KERN_SUCCESS) {
        return 0;
    }
    
    // Free memory = free + inactive + purgeable
    natural_t free_memory = vm_stat.free_count + vm_stat.inactive_count + vm_stat.purgeable_count;
    return free_memory * page_size;
    
#elif defined(_WIN32)
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);
    return memInfo.ullAvailPhys;
    
#elif defined(__linux__)
    struct sysinfo si;
    if (sysinfo(&si) == 0) {
        return si.freeram * si.mem_unit;
    }
    return 0;
#else
    return 0;
#endif
}

CacheManager::CacheManager(size_t initialMaxSize) 
    : mCurrentMaxSize(initialMaxSize),
      mLastResizeCheck(std::chrono::steady_clock::now()) {
    
    mStats.maxSize = initialMaxSize;
    
    // Set default memory callback
    mAvailableMemoryCallback = getAvailableSystemMemory;
}

CacheManager::~CacheManager() {
    clear();
}

CacheManager::AudioBufferPtr CacheManager::get(const AudioBufferKey& key,
                                               std::function<AudioBufferPtr()> loader) {
    std::lock_guard<std::mutex> lock(mCacheMutex);
    
    // Check if item exists in cache
    auto it = mCacheMap.find(key);
    if (it != mCacheMap.end()) {
        // Cache hit
        mStats.hits++;
        
        // Update access info
        auto& node = *(it->second);
        node.lastAccess = std::chrono::steady_clock::now();
        node.accessCount++;
        
        // Move to front of LRU list
        moveToFront(it->second);
        
        mStats.updateHitRate();
        return node.buffer;
    }
    
    // Cache miss
    mStats.misses++;
    mStats.updateHitRate();
    
    // Load the buffer
    AudioBufferPtr buffer = loader();
    if (!buffer) {
        return nullptr;
    }
    
    size_t bufferSize = buffer->size() * sizeof(float);
    
    // Check if we need to evict items
    if (mStats.currentSize + bufferSize > mCurrentMaxSize) {
        evictLRU(bufferSize);
    }
    
    // Add to cache
    CacheNode node;
    node.key = key;
    node.buffer = buffer;
    node.size = bufferSize;
    node.lastAccess = std::chrono::steady_clock::now();
    node.accessCount = 1;
    
    mCacheList.push_front(node);
    mCacheMap[key] = mCacheList.begin();
    
    mStats.currentSize += bufferSize;
    
    // Check if we should adjust cache size
    checkAndAdjustSize();
    
    return buffer;
}

void CacheManager::prefetch(const std::vector<AudioBufferKey>& keys,
                           std::function<AudioBufferPtr(const AudioBufferKey&)> loader) {
    // Launch prefetch in background thread to avoid blocking
    std::thread prefetchThread([this, keys, loader]() {
        for (const auto& key : keys) {
            // Check if already in cache before loading
            {
                std::lock_guard<std::mutex> lock(mCacheMutex);
                if (mCacheMap.find(key) != mCacheMap.end()) {
                    continue;  // Already cached
                }
            }
            
            // Load buffer
            AudioBufferPtr buffer = loader(key);
            if (!buffer) continue;
            
            // Add to cache
            get(key, [buffer]() { return buffer; });
        }
    });
    
    prefetchThread.detach();
}

void CacheManager::clear() {
    std::lock_guard<std::mutex> lock(mCacheMutex);
    
    mCacheList.clear();
    mCacheMap.clear();
    mStats.currentSize = 0;
    mStats.evictions += mCacheMap.size();
}

void CacheManager::adjustCacheSize() {
    std::lock_guard<std::mutex> lock(mCacheMutex);
    checkAndAdjustSize();
}

void CacheManager::setSizeLimits(size_t minSize, size_t maxSize) {
    std::lock_guard<std::mutex> lock(mCacheMutex);
    
    mMinCacheSize = minSize;
    mMaxCacheSize = maxSize;
    
    // Ensure current max is within bounds
    mCurrentMaxSize = std::max(mMinCacheSize, std::min(mCurrentMaxSize, mMaxCacheSize));
    mStats.maxSize = mCurrentMaxSize;
    
    // Evict if necessary
    if (mStats.currentSize > mCurrentMaxSize) {
        evictLRU(mStats.currentSize - mCurrentMaxSize);
    }
}

void CacheManager::evictLRU(size_t requiredSpace) {
    size_t freedSpace = 0;
    
    while (!mCacheList.empty() && freedSpace < requiredSpace) {
        auto& node = mCacheList.back();
        
        freedSpace += node.size;
        mStats.currentSize -= node.size;
        mStats.evictions++;
        
        mCacheMap.erase(node.key);
        mCacheList.pop_back();
    }
}

void CacheManager::moveToFront(CacheList::iterator it) {
    if (it != mCacheList.begin()) {
        mCacheList.splice(mCacheList.begin(), mCacheList, it);
    }
}

void CacheManager::checkAndAdjustSize() {
    auto now = std::chrono::steady_clock::now();
    
    if (now - mLastResizeCheck < mResizeInterval) {
        return;  // Too soon to check again
    }
    
    mLastResizeCheck = now;
    
    size_t optimalSize = calculateOptimalSize();
    
    if (optimalSize != mCurrentMaxSize) {
        mCurrentMaxSize = optimalSize;
        mStats.maxSize = optimalSize;
        
        // Evict if we've shrunk
        if (mStats.currentSize > mCurrentMaxSize) {
            evictLRU(mStats.currentSize - mCurrentMaxSize);
        }
    }
}

size_t CacheManager::calculateOptimalSize() {
    // Get available system memory
    size_t availableMemory = 0;
    if (mAvailableMemoryCallback) {
        availableMemory = mAvailableMemoryCallback();
    }
    
    // Calculate based on hit rate and available memory
    double hitRate = mStats.hitRate.load();
    
    // Base calculation: use up to 25% of available memory
    size_t targetSize = availableMemory / 4;
    
    // Adjust based on hit rate
    if (hitRate > 0.9) {
        // Very high hit rate - increase cache size
        targetSize = static_cast<size_t>(targetSize * 1.2);
    } else if (hitRate > 0.8) {
        // Good hit rate - maintain current size
        targetSize = mCurrentMaxSize;
    } else if (hitRate < 0.5 && mStats.hits + mStats.misses > 100) {
        // Poor hit rate with sufficient samples - decrease size
        targetSize = static_cast<size_t>(targetSize * 0.8);
    }
    
    // Apply bounds
    targetSize = std::max(mMinCacheSize, std::min(targetSize, mMaxCacheSize));
    
    // Don't make small adjustments
    size_t sizeDiff = (targetSize > mCurrentMaxSize) ? 
                      (targetSize - mCurrentMaxSize) : 
                      (mCurrentMaxSize - targetSize);
    
    if (sizeDiff < mCurrentMaxSize * 0.1) {
        return mCurrentMaxSize;  // Less than 10% change - keep current
    }
    
    return targetSize;
}

} // namespace Explorer
} // namespace Acorex