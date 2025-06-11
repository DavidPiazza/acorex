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

#include "AudioTransportN.hpp"
#include <complex>
#include <algorithm>

namespace Acorex {
namespace Explorer {

void AudioTransportN::initializeBuffers(fluid::index maxFFTSize) {
    fluid::index bins = maxFFTSize / 2 + 1;
    
    // Initialize STFT processors
    mSTFT = std::make_unique<fluid::algorithm::STFT>(maxFFTSize, maxFFTSize, maxFFTSize / 2);
    mISTFT = std::make_unique<fluid::algorithm::ISTFT>(maxFFTSize, maxFFTSize, maxFFTSize / 2);
    mReassignSTFT = std::make_unique<fluid::algorithm::STFT>(
        maxFFTSize, maxFFTSize, maxFFTSize / 2,
        static_cast<fluid::index>(fluid::algorithm::WindowFuncs::WindowTypes::kHannD)
    );
    
    // Initialize window
    mWindow = Eigen::ArrayXd(maxFFTSize);
    mWindow.setZero();
    // Create window function
    Eigen::Ref<Eigen::ArrayXd> windowRef(mWindow.head(mWindowSize));
    fluid::algorithm::WindowFuncs::map()[fluid::algorithm::WindowFuncs::WindowTypes::kHann](
        mWindowSize, windowRef
    );
    mWindowSquared = mWindow * mWindow;
    
    // Initialize phase tracking
    mPhase = Eigen::ArrayXd::Zero(bins);
    mBinFreqs = Eigen::ArrayXd::LinSpaced(bins, 0, bins - 1) * (2 * M_PI) / mFFTSize;
    mPhaseDiff = mBinFreqs * mHopSize;
    
    // Pre-allocate buffers for maximum sources
    mSpectraBuffers.reserve(MAX_SOURCES);
    mSpectraDhBuffers.reserve(MAX_SOURCES);
    mMagnitudeBuffers.reserve(MAX_SOURCES);
    mPhaseBuffers.reserve(MAX_SOURCES);
    mReassignedFreqBuffers.reserve(MAX_SOURCES);
    
    for (size_t i = 0; i < MAX_SOURCES; ++i) {
        mSpectraBuffers.emplace_back(bins);
        mSpectraDhBuffers.emplace_back(bins);
        mMagnitudeBuffers.emplace_back(bins);
        mPhaseBuffers.emplace_back(bins);
        mReassignedFreqBuffers.emplace_back(bins);
    }
    
    mAccumulatorSpectrum = Eigen::ArrayXcd(bins);
    mAccumulatorMagnitude = Eigen::ArrayXd(bins);
    mAccumulatorPhase = Eigen::ArrayXd(bins);
}

void AudioTransportN::initN(fluid::index windowSize, fluid::index fftSize, fluid::index hopSize) {
    mWindowSize = windowSize;
    mFFTSize = fftSize;
    mHopSize = hopSize;
    mBins = fftSize / 2 + 1;
    
    // Initialize all transport instances
    for (auto& transport : mTransports) {
        transport->init(windowSize, fftSize, hopSize);
    }
    
    // Resize STFT processors
    mSTFT->resize(windowSize, fftSize, hopSize);
    mISTFT->resize(windowSize, fftSize, hopSize);
    mReassignSTFT->resize(windowSize, fftSize, hopSize);
    
    // Update window
    mWindow.resize(windowSize);
    mWindow.setZero();
    Eigen::Ref<Eigen::ArrayXd> windowRef(mWindow.head(windowSize));
    fluid::algorithm::WindowFuncs::map()[fluid::algorithm::WindowFuncs::WindowTypes::kHann](
        windowSize, windowRef
    );
    mWindowSquared = mWindow * mWindow;
    
    // Update phase tracking
    mPhase = Eigen::ArrayXd::Zero(mBins);
    mBinFreqs = Eigen::ArrayXd::LinSpaced(mBins, 0, mBins - 1) * (2 * M_PI) / fftSize;
    mPhaseDiff = mBinFreqs * hopSize;
    
    // Ensure our buffers are sized correctly for the new FFT size
    for (size_t i = 0; i < MAX_SOURCES; ++i) {
        mSpectraBuffers[i].resize(mBins);
        mSpectraDhBuffers[i].resize(mBins);
        mMagnitudeBuffers[i].resize(mBins);
        mPhaseBuffers[i].resize(mBins);
        mReassignedFreqBuffers[i].resize(mBins);
    }
    
    mAccumulatorSpectrum.resize(mBins);
    mAccumulatorMagnitude.resize(mBins);
    mAccumulatorPhase.resize(mBins);
    
    mInitialized = true;
}

void AudioTransportN::processFrameN(
    const std::vector<fluid::RealVectorView>& frames,
    const BarycentricWeights& weights,
    fluid::RealMatrixView out) {
    
    using namespace fluid::algorithm;
    
    // Validate inputs
    if (frames.size() != weights.size()) {
        throw std::invalid_argument("Number of frames must match number of weights");
    }
    
    if (!weights.isValid()) {
        throw std::invalid_argument("Barycentric weights must sum to 1.0");
    }
    
    if (frames.size() > MAX_SOURCES) {
        throw std::invalid_argument("Too many sources, maximum is " + std::to_string(MAX_SOURCES));
    }
    
    if (!initialized()) {
        throw std::runtime_error("AudioTransportN not initialized");
    }
    
    const size_t nSources = frames.size();
    
    // Handle edge cases
    if (nSources == 0) {
        // Output silence
        for (fluid::index i = 0; i < out.cols(); ++i) {
            out(0, i) = 0.0;
            out(1, i) = i < mWindowSize ? mWindowSquared(i) : 0.0;
        }
        return;
    }
    
    if (nSources == 1) {
        // Single source - just copy through with windowing
        for (fluid::index i = 0; i < out.cols() && i < frames[0].size(); ++i) {
            out(0, i) = frames[0][i];
            out(1, i) = i < mWindowSize ? mWindowSquared(i) : 0.0;
        }
        return;
    }
    
    if (nSources == 2 && mTransports[0]->initialized()) {
        // Use optimized pairwise processing
        mTransports[0]->processFrame(frames[0], frames[1], weights[1], out, mAllocator);
        return;
    }
    
    // N-way processing (N >= 3)
    
    // Process each frame through STFT
    std::vector<double> totalMagnitudes(nSources, 0.0);
    
    for (size_t i = 0; i < nSources; ++i) {
        // Copy frame data to Eigen array
        Eigen::ArrayXd frame(frames[i].size());
        for (fluid::index j = 0; j < frames[i].size(); ++j) {
            frame(j) = frames[i][j];
        }
        
        // Forward transform
        mSTFT->processFrame(frame, mSpectraBuffers[i]);
        mReassignSTFT->processFrame(frame, mSpectraDhBuffers[i]);
        
        // Extract magnitude and phase
        mMagnitudeBuffers[i] = mSpectraBuffers[i].abs();
        mPhaseBuffers[i] = mSpectraBuffers[i].arg();
        totalMagnitudes[i] = mMagnitudeBuffers[i].sum();
        
        // Compute reassigned frequencies
        mReassignedFreqBuffers[i] = mBinFreqs.head(mBins) - 
            (mSpectraDhBuffers[i] / mSpectraBuffers[i]).imag();
    }
    
    // Check if all sources are silent
    bool allSilent = std::all_of(totalMagnitudes.begin(), totalMagnitudes.end(),
                                 [](double m) { return m <= 0.0; });
    
    if (allSilent) {
        // Output silence
        for (fluid::index i = 0; i < out.cols(); ++i) {
            out(0, i) = 0.0;
            out(1, i) = i < mWindowSize ? mWindowSquared(i) : 0.0;
        }
        return;
    }
    
    // Initialize accumulator
    mAccumulatorSpectrum.setZero();
    
    // Compute N-way spectral interpolation
    // For now, use a simplified approach - we'll implement full N-way transport in the next subtask
    
    // Magnitude interpolation using weighted geometric mean
    mAccumulatorMagnitude.setOnes();
    for (size_t i = 0; i < nSources; ++i) {
        if (weights[i] > 0 && totalMagnitudes[i] > 0) {
            // Weighted geometric mean: product(mag_i^weight_i)
            mAccumulatorMagnitude *= (mMagnitudeBuffers[i] + 1e-10).pow(weights[i]);
        }
    }
    
    // Phase interpolation using circular statistics
    Eigen::ArrayXd sinSum = Eigen::ArrayXd::Zero(mBins);
    Eigen::ArrayXd cosSum = Eigen::ArrayXd::Zero(mBins);
    
    for (size_t i = 0; i < nSources; ++i) {
        if (weights[i] > 0 && totalMagnitudes[i] > 0) {
            // Weight phases by both barycentric weight and local magnitude
            Eigen::ArrayXd weightedMag = weights[i] * mMagnitudeBuffers[i];
            sinSum += weightedMag * mPhaseBuffers[i].sin();
            cosSum += weightedMag * mPhaseBuffers[i].cos();
        }
    }
    
    // Compute circular mean phase
    mAccumulatorPhase = sinSum.binaryExpr(cosSum, 
        [](double s, double c) { return std::atan2(s, c); });
    
    // Reconstruct complex spectrum
    for (fluid::index i = 0; i < mBins; ++i) {
        mAccumulatorSpectrum(i) = std::polar(mAccumulatorMagnitude(i), mAccumulatorPhase(i));
    }
    
    // Inverse transform
    Eigen::ArrayXd output(frames[0].size());
    mISTFT->processFrame(mAccumulatorSpectrum, output);
    
    // Write output
    for (fluid::index i = 0; i < out.cols() && i < output.size(); ++i) {
        out(0, i) = output(i);
        out(1, i) = i < mWindowSize ? mWindowSquared(i) : 0.0;
    }
}

void AudioTransportN::computeWeightedGeometricMean(
    const std::vector<Eigen::ArrayXd>& magnitudes,
    const BarycentricWeights& weights,
    Eigen::ArrayXd& result) {
    
    result.setOnes();
    
    for (size_t i = 0; i < magnitudes.size(); ++i) {
        if (weights[i] > 0) {
            // Add small epsilon to avoid log(0)
            result *= (magnitudes[i] + 1e-10).pow(weights[i]);
        }
    }
}

void AudioTransportN::computeWeightedCircularMean(
    const std::vector<Eigen::ArrayXd>& phases,
    const std::vector<Eigen::ArrayXd>& magnitudes,
    const BarycentricWeights& weights,
    Eigen::ArrayXd& result) {
    
    Eigen::ArrayXd sinSum = Eigen::ArrayXd::Zero(result.size());
    Eigen::ArrayXd cosSum = Eigen::ArrayXd::Zero(result.size());
    
    for (size_t i = 0; i < phases.size(); ++i) {
        if (weights[i] > 0) {
            // Weight by both barycentric coordinate and magnitude
            Eigen::ArrayXd weight = weights[i] * magnitudes[i];
            sinSum += weight * phases[i].sin();
            cosSum += weight * phases[i].cos();
        }
    }
    
    // Compute circular mean
    result = sinSum.binaryExpr(cosSum, 
        [](double s, double c) { return std::atan2(s, c); });
}

void AudioTransportN::computeNWayTransport(
    const std::vector<std::vector<SpetralMass>>& masses,
    const BarycentricWeights& weights,
    std::vector<std::tuple<std::vector<fluid::index>, double>>& transport) {
    
    // This is a placeholder for the full N-way optimal transport implementation
    // For now, we'll implement a simplified version that works but isn't optimal
    
    transport.clear();
    
    // TODO: Implement full N-way Wasserstein barycenters
    // This requires solving a more complex optimization problem
    // For now, we use a heuristic approach
    
    if (masses.empty()) return;
    
    // Find the source with maximum weight as reference
    size_t refIdx = 0;
    double maxWeight = weights[0];
    for (size_t i = 1; i < weights.size(); ++i) {
        if (weights[i] > maxWeight) {
            maxWeight = weights[i];
            refIdx = i;
        }
    }
    
    // For each mass in the reference source, find corresponding masses in other sources
    for (size_t m = 0; m < masses[refIdx].size(); ++m) {
        std::vector<fluid::index> indices(masses.size());
        indices[refIdx] = m;
        
        // Simple nearest-neighbor matching for other sources
        for (size_t s = 0; s < masses.size(); ++s) {
            if (s != refIdx && !masses[s].empty()) {
                // Find closest mass by center frequency
                fluid::index refCenter = masses[refIdx][m].centerBin;
                fluid::index bestMatch = 0;
                fluid::index minDist = std::abs(masses[s][0].centerBin - refCenter);
                
                for (size_t i = 1; i < masses[s].size(); ++i) {
                    fluid::index dist = std::abs(masses[s][i].centerBin - refCenter);
                    if (dist < minDist) {
                        minDist = dist;
                        bestMatch = i;
                    }
                }
                
                indices[s] = bestMatch;
            }
        }
        
        // Add transport tuple with combined mass
        double combinedMass = 0.0;
        for (size_t s = 0; s < masses.size(); ++s) {
            if (indices[s] < masses[s].size()) {
                combinedMass += weights[s] * masses[s][indices[s]].mass;
            }
        }
        
        transport.emplace_back(indices, combinedMass);
    }
}

} // namespace Explorer
} // namespace Acorex