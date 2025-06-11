#include "MemoryMappedIO.h"
#include "Data.h"
#include <iostream>
#include <chrono>
#include <random>

// Simple test program for memory-mapped I/O
// This file can be compiled separately for testing

using namespace Utils;

void generateTestData(TransportData& data, size_t fileCount, size_t frameCount, size_t frameSize) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    
    data.clear();
    data.reserveFiles(fileCount);
    
    for (size_t i = 0; i < fileCount; ++i) {
        data.addFile();
        data.reserveFrames(i, frameCount);
        
        for (size_t j = 0; j < frameCount; ++j) {
            TransportFrame frame(frameSize);
            
            // Fill with random data
            for (size_t k = 0; k < frameSize; ++k) {
                frame.magnitude[k] = dis(gen);
                frame.phase[k] = dis(gen) * M_PI;
                frame.dH[k] = dis(gen) * 0.1;
            }
            
            data.frames[i].push_back(std::move(frame));
        }
    }
}

bool compareFrames(const TransportFrame& a, const TransportFrame& b, double tolerance = 1e-10) {
    if (a.frameSize() != b.frameSize()) {
        return false;
    }
    
    for (size_t i = 0; i < a.frameSize(); ++i) {
        if (std::abs(a.magnitude[i] - b.magnitude[i]) > tolerance ||
            std::abs(a.phase[i] - b.phase[i]) > tolerance ||
            std::abs(a.dH[i] - b.dH[i]) > tolerance) {
            return false;
        }
    }
    
    return true;
}

bool testBasicIO() {
    std::cout << "Testing basic memory-mapped I/O..." << std::endl;
    
    // Generate test data
    TransportData originalData;
    generateTestData(originalData, 5, 100, 512);
    
    // Write to file
    TransportDataIO writer;
    if (!writer.write("test_transport.bin", originalData)) {
        std::cerr << "Failed to write test data" << std::endl;
        return false;
    }
    
    // Read back
    TransportData readData;
    TransportDataIO reader;
    if (!reader.read("test_transport.bin", readData)) {
        std::cerr << "Failed to read test data" << std::endl;
        return false;
    }
    
    // Compare
    if (originalData.fileCount() != readData.fileCount()) {
        std::cerr << "File count mismatch" << std::endl;
        return false;
    }
    
    for (size_t i = 0; i < originalData.fileCount(); ++i) {
        if (originalData.frameCount(i) != readData.frameCount(i)) {
            std::cerr << "Frame count mismatch for file " << i << std::endl;
            return false;
        }
        
        for (size_t j = 0; j < originalData.frameCount(i); ++j) {
            const auto* origFrame = originalData.getFrame(i, j);
            const auto* readFrame = readData.getFrame(i, j);
            
            if (!origFrame || !readFrame) {
                std::cerr << "Null frame at [" << i << "][" << j << "]" << std::endl;
                return false;
            }
            
            if (!compareFrames(*origFrame, *readFrame)) {
                std::cerr << "Frame data mismatch at [" << i << "][" << j << "]" << std::endl;
                return false;
            }
        }
    }
    
    std::cout << "Basic I/O test passed!" << std::endl;
    return true;
}

bool testLazyLoading() {
    std::cout << "Testing lazy loading..." << std::endl;
    
    // Use existing file from basic test
    TransportDataIO reader;
    if (!reader.openForReading("test_transport.bin")) {
        std::cerr << "Failed to open file for lazy reading" << std::endl;
        return false;
    }
    
    uint32_t fileCount, totalFrames;
    if (!reader.getMetadata(fileCount, totalFrames)) {
        std::cerr << "Failed to get metadata" << std::endl;
        return false;
    }
    
    std::cout << "Files: " << fileCount << ", Total frames: " << totalFrames << std::endl;
    
    // Load a specific frame
    auto frame = reader.loadFrame(2, 50);
    if (!frame || !frame->isValid()) {
        std::cerr << "Failed to load specific frame" << std::endl;
        return false;
    }
    
    std::cout << "Loaded frame size: " << frame->frameSize() << std::endl;
    std::cout << "Lazy loading test passed!" << std::endl;
    return true;
}

bool testPerformance() {
    std::cout << "Testing performance with large dataset..." << std::endl;
    
    // Generate larger dataset
    TransportData largeData;
    generateTestData(largeData, 50, 1000, 1024);
    
    // Measure write time
    auto start = std::chrono::high_resolution_clock::now();
    TransportDataIO writer;
    if (!writer.write("test_large.bin", largeData)) {
        std::cerr << "Failed to write large dataset" << std::endl;
        return false;
    }
    auto writeTime = std::chrono::high_resolution_clock::now() - start;
    
    // Measure read time
    start = std::chrono::high_resolution_clock::now();
    TransportData readData;
    TransportDataIO reader;
    if (!reader.read("test_large.bin", readData)) {
        std::cerr << "Failed to read large dataset" << std::endl;
        return false;
    }
    auto readTime = std::chrono::high_resolution_clock::now() - start;
    
    auto writeMs = std::chrono::duration_cast<std::chrono::milliseconds>(writeTime).count();
    auto readMs = std::chrono::duration_cast<std::chrono::milliseconds>(readTime).count();
    
    std::cout << "Write time: " << writeMs << " ms" << std::endl;
    std::cout << "Read time: " << readMs << " ms" << std::endl;
    
    // Calculate data size
    size_t totalFrames = 50 * 1000;
    size_t dataSize = totalFrames * 1024 * sizeof(double) * 3;
    double dataSizeMB = dataSize / (1024.0 * 1024.0);
    
    std::cout << "Data size: " << dataSizeMB << " MB" << std::endl;
    std::cout << "Write speed: " << (dataSizeMB / (writeMs / 1000.0)) << " MB/s" << std::endl;
    std::cout << "Read speed: " << (dataSizeMB / (readMs / 1000.0)) << " MB/s" << std::endl;
    
    std::cout << "Performance test passed!" << std::endl;
    return true;
}

#ifdef STANDALONE_TEST
int main() {
    bool allPassed = true;
    
    allPassed &= testBasicIO();
    allPassed &= testLazyLoading();
    allPassed &= testPerformance();
    
    // Clean up test files
    std::remove("test_transport.bin");
    std::remove("test_large.bin");
    
    if (allPassed) {
        std::cout << "\nAll tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << "\nSome tests failed!" << std::endl;
        return 1;
    }
}
#endif