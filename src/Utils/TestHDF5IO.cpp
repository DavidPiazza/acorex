#include "HDF5IO.h"
#include "Data.h"
#include <iostream>
#include <chrono>
#include <random>
#include <filesystem>

// Test program for HDF5 transport data I/O
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

bool testBasicHDF5IO() {
    std::cout << "Testing basic HDF5 I/O..." << std::endl;
    
    // Generate test data
    TransportData originalData;
    generateTestData(originalData, 5, 100, 512);
    
    // Write to HDF5 file with compression
    HDF5TransportIO writer;
    if (!writer.write("test_transport.h5", originalData, 6)) {
        std::cerr << "Failed to write HDF5 test data" << std::endl;
        return false;
    }
    
    std::cout << "Compression ratio: " << writer.getCompressionRatio() << ":1" << std::endl;
    
    // Read back
    TransportData readData;
    HDF5TransportIO reader;
    if (!reader.read("test_transport.h5", readData)) {
        std::cerr << "Failed to read HDF5 test data" << std::endl;
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
    
    std::cout << "Basic HDF5 I/O test passed!" << std::endl;
    return true;
}

bool testHDF5LazyLoading() {
    std::cout << "Testing HDF5 lazy loading..." << std::endl;
    
    // Use existing file from basic test
    HDF5TransportIO reader;
    if (!reader.openForReading("test_transport.h5")) {
        std::cerr << "Failed to open HDF5 file for lazy reading" << std::endl;
        return false;
    }
    
    size_t fileCount, totalFrames;
    if (!reader.getMetadata(fileCount, totalFrames)) {
        std::cerr << "Failed to get HDF5 metadata" << std::endl;
        return false;
    }
    
    std::cout << "Files: " << fileCount << ", Total frames: " << totalFrames << std::endl;
    
    // Load a specific frame
    auto frame = reader.loadFrame(2, 50);
    if (!frame || !frame->isValid()) {
        std::cerr << "Failed to load specific frame from HDF5" << std::endl;
        return false;
    }
    
    std::cout << "Loaded frame size: " << frame->frameSize() << std::endl;
    std::cout << "HDF5 lazy loading test passed!" << std::endl;
    return true;
}

bool testCompressionLevels() {
    std::cout << "Testing HDF5 compression levels..." << std::endl;
    
    // Generate test data
    TransportData testData;
    generateTestData(testData, 10, 200, 1024);
    
    // Calculate uncompressed size
    size_t uncompressedSize = 0;
    for (size_t i = 0; i < testData.fileCount(); ++i) {
        for (size_t j = 0; j < testData.frameCount(i); ++j) {
            uncompressedSize += testData.getFrame(i, j)->frameSize() * sizeof(double) * 3;
        }
    }
    
    std::cout << "Uncompressed size: " << (uncompressedSize / 1024.0 / 1024.0) << " MB" << std::endl;
    
    // Test different compression levels
    for (int level : {0, 1, 3, 6, 9}) {
        HDF5TransportIO writer;
        std::string filename = "test_compress_" + std::to_string(level) + ".h5";
        
        auto start = std::chrono::high_resolution_clock::now();
        if (!writer.write(filename, testData, level)) {
            std::cerr << "Failed to write with compression level " << level << std::endl;
            continue;
        }
        auto writeTime = std::chrono::high_resolution_clock::now() - start;
        
        // Get file size
        size_t fileSize = std::filesystem::file_size(filename);
        double ratio = static_cast<double>(uncompressedSize) / fileSize;
        
        auto writeMs = std::chrono::duration_cast<std::chrono::milliseconds>(writeTime).count();
        
        std::cout << "Level " << level << ": "
                  << "Size=" << (fileSize / 1024.0 / 1024.0) << " MB, "
                  << "Ratio=" << ratio << ":1, "
                  << "Time=" << writeMs << " ms" << std::endl;
        
        // Clean up
        std::filesystem::remove(filename);
    }
    
    std::cout << "Compression levels test passed!" << std::endl;
    return true;
}

bool testPerformanceComparison() {
    std::cout << "Testing performance comparison (Binary vs HDF5)..." << std::endl;
    
    // Generate larger dataset
    TransportData largeData;
    generateTestData(largeData, 50, 1000, 1024);
    
    // Test binary format (from MemoryMappedIO)
    {
        auto start = std::chrono::high_resolution_clock::now();
        TransportDataIO binaryWriter;
        binaryWriter.write("test_binary.bin", largeData);
        auto binaryWriteTime = std::chrono::high_resolution_clock::now() - start;
        
        start = std::chrono::high_resolution_clock::now();
        TransportData binaryReadData;
        TransportDataIO binaryReader;
        binaryReader.read("test_binary.bin", binaryReadData);
        auto binaryReadTime = std::chrono::high_resolution_clock::now() - start;
        
        size_t binarySize = std::filesystem::file_size("test_binary.bin");
        
        auto writeMs = std::chrono::duration_cast<std::chrono::milliseconds>(binaryWriteTime).count();
        auto readMs = std::chrono::duration_cast<std::chrono::milliseconds>(binaryReadTime).count();
        
        std::cout << "Binary format:" << std::endl;
        std::cout << "  Size: " << (binarySize / 1024.0 / 1024.0) << " MB" << std::endl;
        std::cout << "  Write time: " << writeMs << " ms" << std::endl;
        std::cout << "  Read time: " << readMs << " ms" << std::endl;
        
        std::filesystem::remove("test_binary.bin");
    }
    
    // Test HDF5 format
    {
        auto start = std::chrono::high_resolution_clock::now();
        HDF5TransportIO hdf5Writer;
        hdf5Writer.write("test_hdf5.h5", largeData, 6);
        auto hdf5WriteTime = std::chrono::high_resolution_clock::now() - start;
        
        start = std::chrono::high_resolution_clock::now();
        TransportData hdf5ReadData;
        HDF5TransportIO hdf5Reader;
        hdf5Reader.read("test_hdf5.h5", hdf5ReadData);
        auto hdf5ReadTime = std::chrono::high_resolution_clock::now() - start;
        
        size_t hdf5Size = std::filesystem::file_size("test_hdf5.h5");
        
        auto writeMs = std::chrono::duration_cast<std::chrono::milliseconds>(hdf5WriteTime).count();
        auto readMs = std::chrono::duration_cast<std::chrono::milliseconds>(hdf5ReadTime).count();
        
        std::cout << "HDF5 format (compression=6):" << std::endl;
        std::cout << "  Size: " << (hdf5Size / 1024.0 / 1024.0) << " MB" << std::endl;
        std::cout << "  Write time: " << writeMs << " ms" << std::endl;
        std::cout << "  Read time: " << readMs << " ms" << std::endl;
        std::cout << "  Compression ratio: " << hdf5Writer.getCompressionRatio() << ":1" << std::endl;
        
        std::filesystem::remove("test_hdf5.h5");
    }
    
    std::cout << "Performance comparison test passed!" << std::endl;
    return true;
}

#ifdef STANDALONE_TEST
int main() {
    bool allPassed = true;
    
    allPassed &= testBasicHDF5IO();
    allPassed &= testHDF5LazyLoading();
    allPassed &= testCompressionLevels();
    allPassed &= testPerformanceComparison();
    
    // Clean up test files
    std::filesystem::remove("test_transport.h5");
    
    if (allPassed) {
        std::cout << "\nAll HDF5 tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << "\nSome HDF5 tests failed!" << std::endl;
        return 1;
    }
}
#endif