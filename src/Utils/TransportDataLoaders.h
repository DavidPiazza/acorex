#pragma once

#include "TransportFrameCache.h"
#include "HDF5IO.h"
#include "MemoryMappedIO.h"
#include <memory>
#include <string>

namespace Acorex {
namespace Utils {

// HDF5-based TransportData loader
class HDF5TransportDataLoader : public TransportDataLoader {
public:
    explicit HDF5TransportDataLoader(const std::string& filepath) 
        : mIO(std::make_unique<HDF5TransportIO>()) {
        if (!mIO->openForReading(filepath)) {
            throw std::runtime_error("Failed to open HDF5 file: " + filepath);
        }
        mIO->getMetadata(mFileCount, mTotalFrames);
    }
    
    ~HDF5TransportDataLoader() override {
        if (mIO) {
            mIO->close();
        }
    }
    
    std::shared_ptr<TransportFrame> loadFrame(size_t fileIndex, size_t frameIndex) override {
        return mIO->loadFrame(fileIndex, frameIndex);
    }
    
    bool hasFrame(size_t fileIndex, size_t frameIndex) const override {
        // Simple bounds check - could be enhanced with actual metadata
        return fileIndex < mFileCount;
    }
    
private:
    std::unique_ptr<HDF5TransportIO> mIO;
    size_t mFileCount = 0;
    size_t mTotalFrames = 0;
};

// Memory-mapped TransportData loader
class MemoryMappedTransportDataLoader : public TransportDataLoader {
public:
    explicit MemoryMappedTransportDataLoader(const std::string& filepath)
        : mIO(std::make_unique<TransportDataIO>()) {
        if (!mIO->openForReading(filepath)) {
            throw std::runtime_error("Failed to open memory-mapped file: " + filepath);
        }
        mIO->getMetadata(mFileCount, mTotalFrames);
    }
    
    ~MemoryMappedTransportDataLoader() override {
        if (mIO) {
            mIO->close();
        }
    }
    
    std::shared_ptr<TransportFrame> loadFrame(size_t fileIndex, size_t frameIndex) override {
        return mIO->loadFrame(fileIndex, frameIndex);
    }
    
    bool hasFrame(size_t fileIndex, size_t frameIndex) const override {
        return fileIndex < mFileCount;
    }
    
private:
    std::unique_ptr<TransportDataIO> mIO;
    uint32_t mFileCount = 0;
    uint32_t mTotalFrames = 0;
};

// In-memory TransportData loader (for testing or small datasets)
class InMemoryTransportDataLoader : public TransportDataLoader {
public:
    explicit InMemoryTransportDataLoader(std::shared_ptr<TransportData> data)
        : mData(data) {
        if (!mData) {
            throw std::invalid_argument("TransportData cannot be null");
        }
    }
    
    std::shared_ptr<TransportFrame> loadFrame(size_t fileIndex, size_t frameIndex) override {
        auto* frame = mData->getFrame(fileIndex, frameIndex);
        if (frame) {
            // Create a copy to return as shared_ptr
            return std::make_shared<TransportFrame>(*frame);
        }
        return nullptr;
    }
    
    bool hasFrame(size_t fileIndex, size_t frameIndex) const override {
        return mData->hasFile(fileIndex) && frameIndex < mData->frameCount(fileIndex);
    }
    
private:
    std::shared_ptr<TransportData> mData;
};

// Composite loader that tries multiple sources
class CompositeTransportDataLoader : public TransportDataLoader {
public:
    void addLoader(std::shared_ptr<TransportDataLoader> loader) {
        mLoaders.push_back(loader);
    }
    
    std::shared_ptr<TransportFrame> loadFrame(size_t fileIndex, size_t frameIndex) override {
        for (auto& loader : mLoaders) {
            if (loader->hasFrame(fileIndex, frameIndex)) {
                auto frame = loader->loadFrame(fileIndex, frameIndex);
                if (frame) {
                    return frame;
                }
            }
        }
        return nullptr;
    }
    
    bool hasFrame(size_t fileIndex, size_t frameIndex) const override {
        return std::any_of(mLoaders.begin(), mLoaders.end(),
            [fileIndex, frameIndex](const auto& loader) {
                return loader->hasFrame(fileIndex, frameIndex);
            });
    }
    
private:
    std::vector<std::shared_ptr<TransportDataLoader>> mLoaders;
};

// Factory function to create appropriate loader based on file extension
inline std::shared_ptr<TransportDataLoader> createTransportDataLoader(const std::string& filepath) {
    if (filepath.find(".h5") != std::string::npos || 
        filepath.find(".hdf5") != std::string::npos) {
        return std::make_shared<HDF5TransportDataLoader>(filepath);
    } else if (filepath.find(".bin") != std::string::npos ||
               filepath.find(".dat") != std::string::npos) {
        return std::make_shared<MemoryMappedTransportDataLoader>(filepath);
    } else {
        throw std::runtime_error("Unsupported file format: " + filepath);
    }
}

} // namespace Utils
} // namespace Acorex