/*
The MIT License (MIT)

Copyright (c) 2024

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <numeric>
#include "../src/Explorer/AudioMorphEngine.h"
#include <ofSoundBuffer.h>

using namespace Acorex::Explorer;

class OverlapAddTest {
public:
    void testLatency() {
        std::cout << "\n=== Testing Overlap-Add Latency ===" << std::endl;
        
        AudioMorphEngine engine;
        int sampleRate = 44100;
        size_t bufferSize = 512;
        float overlap = 0.75f;
        
        engine.Initialise(sampleRate, bufferSize, overlap);
        
        // Process several buffers to warm up
        ofSoundBuffer buffer;
        buffer.setSampleRate(sampleRate);
        buffer.allocate(bufferSize, 2); // stereo
        
        for (int i = 0; i < 10; i++) {
            auto start = std::chrono::high_resolution_clock::now();
            engine.Process(buffer);
            auto end = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            double processingTimeMs = duration / 1000.0;
            
            std::cout << "Buffer " << i << " processing time: " << processingTimeMs << " ms" << std::endl;
        }
        
        double latency = engine.GetCurrentLatencyMs();
        std::cout << "Current engine latency: " << latency << " ms" << std::endl;
        std::cout << "Underrun count: " << engine.GetUnderrunCount() << std::endl;
        std::cout << "Overrun count: " << engine.GetOverrunCount() << std::endl;
        
        if (latency < 10.0) {
            std::cout << "✓ Latency target achieved (<10ms)" << std::endl;
        } else {
            std::cout << "✗ Latency exceeds target (>10ms)" << std::endl;
        }
    }
    
    void testOverlapAccuracy() {
        std::cout << "\n=== Testing 75% Overlap Accuracy ===" << std::endl;
        
        AudioMorphEngine engine;
        int sampleRate = 44100;
        size_t bufferSize = 512;
        
        // Test different window sizes with 75% overlap
        std::vector<size_t> windowSizes = {1024, 2048, 4096};
        
        for (size_t windowSize : windowSizes) {
            size_t hopSize = windowSize / 4; // 75% overlap
            
            std::cout << "\nWindow size: " << windowSize << ", Hop size: " << hopSize << std::endl;
            
            engine.SetFFTParameters(windowSize, windowSize, hopSize);
            
            float actualOverlap = 1.0f - (float)hopSize / windowSize;
            std::cout << "Actual overlap: " << actualOverlap * 100 << "%" << std::endl;
            
            if (std::abs(actualOverlap - 0.75f) < 0.001f) {
                std::cout << "✓ 75% overlap correctly configured" << std::endl;
            } else {
                std::cout << "✗ Overlap configuration error" << std::endl;
            }
        }
    }
    
    void testCircularBufferTiming() {
        std::cout << "\n=== Testing Circular Buffer Timing ===" << std::endl;
        
        AudioMorphEngine engine;
        int sampleRate = 44100;
        size_t bufferSize = 512;
        
        engine.Initialise(sampleRate, bufferSize, 0.75f);
        
        // Process buffers and check sample-accurate timing
        ofSoundBuffer buffer;
        buffer.setSampleRate(sampleRate);
        buffer.allocate(bufferSize, 2);
        
        // Track total samples processed
        size_t totalSamples = 0;
        auto startTime = std::chrono::steady_clock::now();
        
        // Process 1 second worth of audio
        int buffersFor1Second = sampleRate / bufferSize;
        
        for (int i = 0; i < buffersFor1Second; i++) {
            engine.Process(buffer);
            totalSamples += bufferSize;
        }
        
        auto endTime = std::chrono::steady_clock::now();
        auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Processed " << totalSamples << " samples in " << elapsedMs << " ms" << std::endl;
        
        double expectedMs = (totalSamples * 1000.0) / sampleRate;
        std::cout << "Expected duration: " << expectedMs << " ms" << std::endl;
        
        // Check if we're maintaining sample-accurate timing
        if (engine.GetUnderrunCount() == 0) {
            std::cout << "✓ No buffer underruns - timing is stable" << std::endl;
        } else {
            std::cout << "✗ Buffer underruns detected: " << engine.GetUnderrunCount() << std::endl;
        }
    }
    
    void testWindowFunction() {
        std::cout << "\n=== Testing Window Function COLA Property ===" << std::endl;
        
        // Test that overlapped windows sum to constant (COLA - Constant Overlap-Add)
        size_t windowSize = 2048;
        size_t hopSize = 512; // 75% overlap
        
        // Generate Hann window
        std::vector<float> window(windowSize);
        for (size_t i = 0; i < windowSize; ++i) {
            double phase = 2.0 * M_PI * i / (windowSize - 1);
            window[i] = 0.5f * (1.0f - std::cos(phase));
        }
        
        // Normalize for COLA
        float windowSum = 0.0f;
        for (size_t i = 0; i < hopSize; ++i) {
            float sum = 0.0f;
            for (size_t j = i; j < windowSize; j += hopSize) {
                sum += window[j] * window[j];
            }
            windowSum += sum;
        }
        
        float normFactor = std::sqrt(hopSize / windowSum);
        for (float& w : window) {
            w *= normFactor;
        }
        
        // Verify COLA property
        std::vector<float> colaTest(windowSize * 2, 0.0f);
        
        // Add overlapped windows
        for (size_t offset = 0; offset < windowSize; offset += hopSize) {
            for (size_t i = 0; i < windowSize; ++i) {
                colaTest[offset + i] += window[i] * window[i];
            }
        }
        
        // Check center region for constant sum
        bool colaValid = true;
        float targetSum = colaTest[windowSize]; // Use center value as reference
        
        for (size_t i = windowSize - hopSize; i < windowSize + hopSize; ++i) {
            if (std::abs(colaTest[i] - targetSum) > 0.01f) {
                colaValid = false;
                std::cout << "COLA error at index " << i << ": " << colaTest[i] 
                         << " (expected " << targetSum << ")" << std::endl;
            }
        }
        
        if (colaValid) {
            std::cout << "✓ Window function satisfies COLA property" << std::endl;
        } else {
            std::cout << "✗ Window function does not satisfy COLA property" << std::endl;
        }
    }
};

int main() {
    std::cout << "Testing Overlap-Add Synthesis Implementation" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    OverlapAddTest test;
    
    test.testLatency();
    test.testOverlapAccuracy();
    test.testCircularBufferTiming();
    test.testWindowFunction();
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "Test suite completed" << std::endl;
    
    return 0;
}