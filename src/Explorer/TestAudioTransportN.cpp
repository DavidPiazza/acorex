/*
Test program for AudioTransportN N-way optimal transport functionality
*/

#include "AudioTransportN.hpp"
#include <iostream>
#include <vector>
#include <random>

using namespace Acorex::Explorer;

void testNWayTransport() {
    // Create a simple allocator
    fluid::RealTimeMemory allocator(1024 * 1024); // 1MB
    
    // Create AudioTransportN instance
    const fluid::index maxFFTSize = 2048;
    AudioTransportN transport(maxFFTSize, allocator);
    
    // Initialize with typical values
    const fluid::index windowSize = 1024;
    const fluid::index fftSize = 1024;
    const fluid::index hopSize = 512;
    transport.initN(windowSize, fftSize, hopSize);
    
    // Create test audio frames (3 sources)
    const int frameSize = 512;
    const int nSources = 3;
    
    std::vector<std::vector<double>> testFrames(nSources);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.5, 0.5);
    
    // Generate test signals
    for (int s = 0; s < nSources; ++s) {
        testFrames[s].resize(frameSize);
        // Create simple sine waves at different frequencies
        double freq = 440.0 * (s + 1); // 440Hz, 880Hz, 1320Hz
        double phase = 0.0;
        double phaseInc = 2.0 * M_PI * freq / 44100.0;
        
        for (int i = 0; i < frameSize; ++i) {
            testFrames[s][i] = 0.5 * std::sin(phase) + 0.1 * dis(gen);
            phase += phaseInc;
        }
    }
    
    // Create fluid views
    std::vector<fluid::RealVectorView> frameViews;
    for (const auto& frame : testFrames) {
        frameViews.emplace_back(const_cast<double*>(frame.data()), 0, frameSize);
    }
    
    // Test different weight configurations
    std::vector<AudioTransportN::BarycentricWeights> testWeights = {
        {1.0, 0.0, 0.0},  // 100% source 1
        {0.0, 1.0, 0.0},  // 100% source 2
        {0.0, 0.0, 1.0},  // 100% source 3
        {0.5, 0.5, 0.0},  // 50/50 sources 1&2
        {0.33, 0.33, 0.34}, // Equal mix
        {0.7, 0.2, 0.1}   // Weighted mix
    };
    
    // Output buffer
    std::vector<double> outputData(frameSize * 2); // audio + window
    fluid::RealMatrixView output(outputData.data(), 0, 2, frameSize);
    
    std::cout << "Testing N-way AudioTransport with " << nSources << " sources\n";
    std::cout << "Frame size: " << frameSize << ", FFT size: " << fftSize << "\n\n";
    
    for (size_t w = 0; w < testWeights.size(); ++w) {
        std::cout << "Test " << (w+1) << " - Weights: [";
        for (size_t i = 0; i < testWeights[w].size(); ++i) {
            std::cout << testWeights[w][i];
            if (i < testWeights[w].size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
        
        try {
            transport.processFrameN(frameViews, testWeights[w], output);
            
            // Calculate output energy
            double energy = 0.0;
            for (int i = 0; i < frameSize; ++i) {
                energy += output(0, i) * output(0, i);
            }
            energy = std::sqrt(energy / frameSize);
            
            std::cout << "  Output RMS energy: " << energy << "\n";
            std::cout << "  Success!\n\n";
            
        } catch (const std::exception& e) {
            std::cout << "  Error: " << e.what() << "\n\n";
        }
    }
    
    // Test edge cases
    std::cout << "Testing edge cases:\n";
    
    // Test with 2 sources (should use optimized path)
    std::vector<fluid::RealVectorView> twoFrames(frameViews.begin(), frameViews.begin() + 2);
    AudioTransportN::BarycentricWeights twoWeights{0.3, 0.7};
    
    std::cout << "2-source test (should use optimized pairwise path): ";
    try {
        transport.processFrameN(twoFrames, twoWeights, output);
        std::cout << "Success\n";
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
    }
    
    // Test with invalid weights
    std::cout << "Invalid weights test (should throw): ";
    AudioTransportN::BarycentricWeights invalidWeights{0.5, 0.3, 0.1}; // Sum != 1
    try {
        transport.processFrameN(frameViews, invalidWeights, output);
        std::cout << "Failed - did not throw!\n";
    } catch (const std::exception& e) {
        std::cout << "Success - caught: " << e.what() << "\n";
    }
    
    std::cout << "\nAll tests completed.\n";
}

int main() {
    try {
        testNWayTransport();
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}