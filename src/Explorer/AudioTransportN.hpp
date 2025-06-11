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

#include <algorithms/public/AudioTransport.hpp>
#include <algorithms/public/STFT.hpp>
#include <algorithms/public/WindowFuncs.hpp>
#include <data/TensorTypes.hpp>
#include <data/FluidMemory.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include "SIMDUtils.hpp"

namespace Acorex {
namespace Explorer {

/**
 * N-way barycentric coordinate system for audio morphing
 * Extends FluCoMa's AudioTransport to support multi-source interpolation
 */
class AudioTransportN {
public:
    using SpetralMass = fluid::algorithm::SpetralMass;
    static constexpr size_t MAX_SOURCES = 8;  // Maximum number of simultaneous sources
    
    /**
     * Barycentric weight structure
     * Weights must sum to 1.0 for valid barycentric coordinates
     */
    struct BarycentricWeights {
        std::vector<double> weights;
        
        BarycentricWeights() = default;
        
        explicit BarycentricWeights(size_t n) : weights(n, 0.0) {
            if (n > 0) weights[0] = 1.0;  // Default to first source
        }
        
        BarycentricWeights(std::initializer_list<double> init) : weights(init) {
            normalize();
        }
        
        void normalize() {
            double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
            if (sum > 0) {
                for (auto& w : weights) w /= sum;
            } else if (!weights.empty()) {
                weights[0] = 1.0;  // Fallback to first source
            }
        }
        
        bool isValid() const {
            if (weights.empty()) return false;
            double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
            return std::abs(sum - 1.0) < 1e-6;
        }
        
        size_t size() const { return weights.size(); }
        
        double& operator[](size_t i) { return weights[i]; }
        const double& operator[](size_t i) const { return weights[i]; }
        
        // Convenience methods for std::vector-like usage
        void push_back(double w) { weights.push_back(w); }
        auto begin() { return weights.begin(); }
        auto end() { return weights.end(); }
        auto begin() const { return weights.begin(); }
        auto end() const { return weights.end(); }
    };
    
    /**
     * Constructor
     */
    AudioTransportN(fluid::index maxFFTSize, fluid::Allocator& alloc)
        : mAllocator(alloc) {
        // Create transport instances for pairwise morphing
        for (size_t i = 0; i < MAX_SOURCES; ++i) {
            mTransports.emplace_back(
                std::make_unique<fluid::algorithm::AudioTransport>(maxFFTSize, alloc)
            );
        }
        initializeBuffers(maxFFTSize);
    }
    
    /**
     * Process N-way interpolation
     * @param frames Array of input frames (must match weights.size())
     * @param weights Barycentric weights (must sum to 1.0)
     * @param out Output matrix [audio, window]
     */
    void processFrameN(const std::vector<fluid::RealVectorView>& frames,
                      const BarycentricWeights& weights,
                      fluid::RealMatrixView out);
    
    /**
     * Process N-way interpolation using weighted geometric mean
     * Simpler alternative to transport-based morphing
     * @param frames Array of input frames (must match weights.size())
     * @param weights Barycentric weights (must sum to 1.0)
     * @param out Output matrix [audio, window]
     */
    void processFrameNGeometric(const std::vector<fluid::RealVectorView>& frames,
                               const BarycentricWeights& weights,
                               fluid::RealMatrixView out);
    
    /**
     * Process N-way interpolation with pre-allocated frame views
     * More efficient for real-time processing
     */
    template<size_t N>
    void processFrameN(const std::array<fluid::RealVectorView, N>& frames,
                      const std::array<double, N>& weights,
                      fluid::RealMatrixView out);
    
    /**
     * Initialize for N-way processing
     * Extends base init() to allocate additional buffers
     */
    void initN(fluid::index windowSize, fluid::index fftSize, fluid::index hopSize);
    
    /**
     * Get maximum supported sources
     */
    static constexpr size_t getMaxSources() { return MAX_SOURCES; }
    
    /**
     * Check if initialized
     */
    bool initialized() const { return mInitialized; }
    
private:
    fluid::Allocator& mAllocator;
    fluid::index mWindowSize{1024};
    fluid::index mFFTSize{1024};
    fluid::index mHopSize{512};
    fluid::index mBins{513};
    bool mInitialized{false};
    
    // Composition: AudioTransport instances for pairwise operations
    std::vector<std::unique_ptr<fluid::algorithm::AudioTransport>> mTransports;
    
    // STFT processors for N-way analysis
    std::unique_ptr<fluid::algorithm::STFT> mSTFT;
    std::unique_ptr<fluid::algorithm::ISTFT> mISTFT;
    std::unique_ptr<fluid::algorithm::STFT> mReassignSTFT;
    
    // Window function
    Eigen::ArrayXd mWindow;
    Eigen::ArrayXd mWindowSquared;
    
    // Intermediate buffers for N-way processing
    std::vector<Eigen::ArrayXcd> mSpectraBuffers;
    std::vector<Eigen::ArrayXcd> mSpectraDhBuffers;
    std::vector<Eigen::ArrayXd> mMagnitudeBuffers;
    std::vector<Eigen::ArrayXd> mPhaseBuffers;
    std::vector<Eigen::ArrayXd> mReassignedFreqBuffers;
    
    // Temporary computation buffers
    Eigen::ArrayXcd mAccumulatorSpectrum;
    Eigen::ArrayXd mAccumulatorMagnitude;
    Eigen::ArrayXd mAccumulatorPhase;
    
    // Phase tracking
    Eigen::ArrayXd mPhase;
    Eigen::ArrayXd mPhaseDiff;
    Eigen::ArrayXd mBinFreqs;
    
    void initializeBuffers(fluid::index maxFFTSize);
    
    /**
     * Compute weighted geometric mean of magnitudes
     * Used for perceptually smooth magnitude interpolation
     */
    void computeWeightedGeometricMean(
        const std::vector<Eigen::ArrayXd>& magnitudes,
        const BarycentricWeights& weights,
        Eigen::ArrayXd& result);
    
    /**
     * Compute weighted circular mean of phases
     * Uses von Mises distribution and circular statistics
     * Handles phase wrapping and discontinuities correctly
     */
    void computeWeightedCircularMean(
        const std::vector<Eigen::ArrayXd>& phases,
        const std::vector<Eigen::ArrayXd>& magnitudes,
        const BarycentricWeights& weights,
        Eigen::ArrayXd& result);
    
    /**
     * Apply phase smoothing to reduce discontinuities
     * Uses local phase coherence to adaptively smooth phase values
     */
    void applyPhaseSmoothing(Eigen::ArrayXd& phases);
    
    /**
     * Compute phase coherence measure between sources
     * Returns value between 0 (incoherent) and 1 (coherent)
     */
    double computePhaseCoherence(
        const std::vector<Eigen::ArrayXd>& phases,
        const std::vector<Eigen::ArrayXd>& magnitudes,
        const BarycentricWeights& weights,
        fluid::index bin);
    
    /**
     * Extended transport matrix computation for N sources
     * Generalizes the pairwise optimal transport to N-way
     */
    void computeNWayTransport(
        const std::vector<std::vector<SpetralMass>>& masses,
        const BarycentricWeights& weights,
        std::vector<std::tuple<std::vector<fluid::index>, double>>& transport);
    
    /**
     * Segment spectrum into spectral masses
     * Adapted from FluCoMa's AudioTransport
     */
    std::vector<SpetralMass> segmentSpectrum(
        const Eigen::Ref<Eigen::ArrayXd> mag,
        const Eigen::Ref<Eigen::ArrayXd> reassignedFreq,
        fluid::Allocator& alloc);
    
    /**
     * Place a spectral mass at a given bin location
     * Adapted from FluCoMa's AudioTransport
     */
    void placeMass(const SpetralMass mass, fluid::index bin, double scale,
                   double centerPhase, const Eigen::Ref<Eigen::ArrayXcd> input, 
                   Eigen::Ref<Eigen::ArrayXcd> output,
                   double nextPhase, Eigen::Ref<Eigen::ArrayXd> amplitudes, 
                   Eigen::Ref<Eigen::ArrayXd> phases);
};

// Template implementation for fixed-size processing
template<size_t N>
void AudioTransportN::processFrameN(
    const std::array<fluid::RealVectorView, N>& frames,
    const std::array<double, N>& weights,
    fluid::RealMatrixView out) {
    
    // Convert to dynamic structures
    std::vector<fluid::RealVectorView> frameVec(frames.begin(), frames.end());
    BarycentricWeights weightVec;
    weightVec.weights.assign(weights.begin(), weights.end());
    
    processFrameN(frameVec, weightVec, out);
}

} // namespace Explorer
} // namespace Acorex