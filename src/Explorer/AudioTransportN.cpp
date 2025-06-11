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
    std::vector<std::vector<SpetralMass>> spectralMasses(nSources);
    
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
        
        // Segment spectrum into masses
        if (totalMagnitudes[i] > 0) {
            spectralMasses[i] = segmentSpectrum(mMagnitudeBuffers[i], 
                                               mReassignedFreqBuffers[i], 
                                               mAllocator);
        }
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
    
    // Compute N-way optimal transport
    std::vector<std::tuple<std::vector<fluid::index>, double>> transportPlan;
    computeNWayTransport(spectralMasses, weights, transportPlan);
    
    // Pre-allocate phase tracking arrays
    Eigen::ArrayXd newAmplitudes = Eigen::ArrayXd::Zero(mBins);
    Eigen::ArrayXd newPhases = Eigen::ArrayXd::Zero(mBins);
    
    // Apply transport plan
    for (const auto& transport : transportPlan) {
        const auto& indices = std::get<0>(transport);
        const double transportMass = std::get<1>(transport);
        
        // Compute barycentric position for this transport
        double interpolatedBin = 0.0;
        double totalWeight = 0.0;
        
        for (size_t s = 0; s < nSources; ++s) {
            if (indices[s] >= 0 && indices[s] < spectralMasses[s].size()) {
                interpolatedBin += weights[s] * spectralMasses[s][indices[s]].centerBin;
                totalWeight += weights[s];
            }
        }
        
        if (totalWeight > 0) {
            interpolatedBin /= totalWeight;
        }
        
        fluid::index targetBin = std::lrint(interpolatedBin);
        if (targetBin < 0) targetBin = 0;
        if (targetBin >= mBins) targetBin = mBins - 1;
        
        // Compute interpolated reassigned frequency
        double interpolatedFreq = 0.0;
        for (size_t s = 0; s < nSources; ++s) {
            if (indices[s] >= 0 && indices[s] < spectralMasses[s].size() && weights[s] > 0) {
                fluid::index centerBin = spectralMasses[s][indices[s]].centerBin;
                if (centerBin >= 0 && centerBin < mBins) {
                    interpolatedFreq += weights[s] * mReassignedFreqBuffers[s](centerBin);
                }
            }
        }
        
        // Phase tracking
        double nextPhase = mPhase(targetBin) + interpolatedFreq * mHopSize;
        double centerPhase = nextPhase - mPhaseDiff(targetBin);
        
        // Place masses from all sources
        for (size_t s = 0; s < nSources; ++s) {
            if (indices[s] >= 0 && indices[s] < spectralMasses[s].size() && weights[s] > 0) {
                const SpetralMass& mass = spectralMasses[s][indices[s]];
                double scale = weights[s] * transportMass / mass.mass;
                
                placeMass(mass, targetBin, scale, centerPhase,
                         mSpectraBuffers[s], mAccumulatorSpectrum, 
                         nextPhase, newAmplitudes, newPhases);
            }
        }
    }
    
    // Update phase tracking
    mPhase = newPhases;
    
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
    
    // Prepare data for SIMD processing
    std::vector<const double*> magPtrs;
    magPtrs.reserve(magnitudes.size());
    for (const auto& mag : magnitudes) {
        magPtrs.push_back(mag.data());
    }
    
    // Use SIMD-optimized implementation
    SIMD::computeWeightedGeometricMeanSIMD(
        magPtrs.data(),
        weights.weights.data(),
        magnitudes.size(),
        result.size(),
        result.data()
    );
}

void AudioTransportN::processFrameNGeometric(
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
    
    // N-way geometric mean processing (N >= 3)
    
    // Process each frame through STFT
    for (size_t i = 0; i < nSources; ++i) {
        // Copy frame data to Eigen array
        Eigen::ArrayXd frame(frames[i].size());
        for (fluid::index j = 0; j < frames[i].size(); ++j) {
            frame(j) = frames[i][j];
        }
        
        // Forward transform
        mSTFT->processFrame(frame, mSpectraBuffers[i]);
        
        // Extract magnitude and phase
        mMagnitudeBuffers[i] = mSpectraBuffers[i].abs();
        mPhaseBuffers[i] = mSpectraBuffers[i].arg();
    }
    
    // Compute weighted geometric mean of magnitudes
    Eigen::ArrayXd interpolatedMagnitude(mBins);
    computeWeightedGeometricMean(mMagnitudeBuffers, weights, interpolatedMagnitude);
    
    // Compute weighted circular mean of phases
    Eigen::ArrayXd interpolatedPhase(mBins);
    computeWeightedCircularMean(mPhaseBuffers, mMagnitudeBuffers, weights, interpolatedPhase);
    
    // Reconstruct spectrum from interpolated magnitude and phase
    mAccumulatorSpectrum = Eigen::ArrayXcd::Zero(mBins);
    for (fluid::index i = 0; i < mBins; ++i) {
        mAccumulatorSpectrum(i) = std::polar(interpolatedMagnitude(i), interpolatedPhase(i));
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

void AudioTransportN::computeWeightedCircularMean(
    const std::vector<Eigen::ArrayXd>& phases,
    const std::vector<Eigen::ArrayXd>& magnitudes,
    const BarycentricWeights& weights,
    Eigen::ArrayXd& result) {
    
    // Prepare data for SIMD processing
    std::vector<const double*> phasePtrs;
    std::vector<const double*> magPtrs;
    phasePtrs.reserve(phases.size());
    magPtrs.reserve(magnitudes.size());
    
    for (const auto& phase : phases) {
        phasePtrs.push_back(phase.data());
    }
    for (const auto& mag : magnitudes) {
        magPtrs.push_back(mag.data());
    }
    
    // Use SIMD-optimized implementation
    SIMD::computeWeightedCircularMeanSIMD(
        phasePtrs.data(),
        magPtrs.data(),
        weights.weights.data(),
        phases.size(),
        result.size(),
        result.data()
    );
    
    // Post-process for phase continuity
    applyPhaseSmoothing(result);
}

void AudioTransportN::applyPhaseSmoothing(Eigen::ArrayXd& phases) {
    // Use SIMD-optimized phase smoothing
    SIMD::applyPhaseSmoothingSIMD(phases.data(), phases.size());
}

double AudioTransportN::computePhaseCoherence(
    const std::vector<Eigen::ArrayXd>& phases,
    const std::vector<Eigen::ArrayXd>& magnitudes,
    const BarycentricWeights& weights,
    fluid::index bin) {
    
    // Compute phase coherence using circular variance
    // Returns 1 for perfectly coherent phases, 0 for random phases
    
    const double epsilon = 1e-10;
    double sumX = 0.0;
    double sumY = 0.0;
    double totalWeight = 0.0;
    
    for (size_t i = 0; i < phases.size(); ++i) {
        if (weights[i] > epsilon && magnitudes[i](bin) > epsilon) {
            double w = weights[i] * magnitudes[i](bin);
            double phase = phases[i](bin);
            
            sumX += w * std::cos(phase);
            sumY += w * std::sin(phase);
            totalWeight += w;
        }
    }
    
    if (totalWeight > epsilon) {
        // Mean resultant length (R)
        double R = std::sqrt(sumX * sumX + sumY * sumY) / totalWeight;
        return R;  // R ranges from 0 (incoherent) to 1 (coherent)
    }
    
    return 0.0;  // No data - assume incoherent
}

void AudioTransportN::computeNWayTransport(
    const std::vector<std::vector<SpetralMass>>& masses,
    const BarycentricWeights& weights,
    std::vector<std::tuple<std::vector<fluid::index>, double>>& transport) {
    
    transport.clear();
    
    if (masses.empty() || masses[0].empty()) return;
    
    const size_t nSources = masses.size();
    
    // Build unified mass list with source indices
    struct IndexedMass {
        fluid::index sourceIdx;
        fluid::index massIdx;
        fluid::index centerBin;
        double mass;
    };
    
    std::vector<IndexedMass> allMasses;
    for (size_t s = 0; s < nSources; ++s) {
        for (size_t m = 0; m < masses[s].size(); ++m) {
            allMasses.push_back({
                static_cast<fluid::index>(s),
                static_cast<fluid::index>(m),
                masses[s][m].centerBin,
                masses[s][m].mass * weights[s]  // Pre-weight the masses
            });
        }
    }
    
    // Sort by center frequency for efficient clustering
    std::sort(allMasses.begin(), allMasses.end(), 
        [](const IndexedMass& a, const IndexedMass& b) {
            return a.centerBin < b.centerBin;
        });
    
    // Compute barycentric clusters using iterative refinement
    // This implements a simplified version of the Wasserstein barycenter algorithm
    const int MAX_ITERATIONS = 10;
    const double CONVERGENCE_THRESHOLD = 1e-6;
    
    // Initialize cluster centers based on weighted average of all masses
    std::vector<double> clusterCenters;
    std::vector<std::vector<size_t>> clusterMembers;
    
    // Use adaptive clustering based on frequency distribution
    
    // Create initial clusters by partitioning the frequency space
    const int nClusters = std::min(static_cast<int>(allMasses.size()), 
                                  static_cast<int>(masses[0].size() * 2));
    
    clusterCenters.resize(nClusters);
    clusterMembers.resize(nClusters);
    
    // Initialize cluster centers uniformly across frequency range
    if (nClusters > 1) {
        fluid::index minBin = allMasses.front().centerBin;
        fluid::index maxBin = allMasses.back().centerBin;
        for (int c = 0; c < nClusters; ++c) {
            clusterCenters[c] = minBin + (maxBin - minBin) * c / (nClusters - 1.0);
        }
    } else {
        clusterCenters[0] = allMasses[allMasses.size() / 2].centerBin;
    }
    
    // Iterative refinement using Lloyd's algorithm adapted for Wasserstein space
    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        // Clear cluster assignments
        for (auto& members : clusterMembers) {
            members.clear();
        }
        
        // Assign masses to nearest cluster center
        for (size_t i = 0; i < allMasses.size(); ++i) {
            int bestCluster = 0;
            double minDist = std::abs(allMasses[i].centerBin - clusterCenters[0]);
            
            for (int c = 1; c < nClusters; ++c) {
                double dist = std::abs(allMasses[i].centerBin - clusterCenters[c]);
                if (dist < minDist) {
                    minDist = dist;
                    bestCluster = c;
                }
            }
            
            clusterMembers[bestCluster].push_back(i);
        }
        
        // Update cluster centers as weighted average
        double maxShift = 0.0;
        for (int c = 0; c < nClusters; ++c) {
            if (clusterMembers[c].empty()) continue;
            
            double weightedSum = 0.0;
            double totalMass = 0.0;
            
            for (size_t idx : clusterMembers[c]) {
                weightedSum += allMasses[idx].centerBin * allMasses[idx].mass;
                totalMass += allMasses[idx].mass;
            }
            
            if (totalMass > 0) {
                double newCenter = weightedSum / totalMass;
                maxShift = std::max(maxShift, std::abs(newCenter - clusterCenters[c]));
                clusterCenters[c] = newCenter;
            }
        }
        
        // Check convergence
        if (maxShift < CONVERGENCE_THRESHOLD) {
            break;
        }
    }
    
    // Build transport plan from clusters
    for (int c = 0; c < nClusters; ++c) {
        if (clusterMembers[c].empty()) continue;
        
        // For each cluster, create a transport tuple
        std::vector<fluid::index> sourceIndices(nSources, -1);
        double clusterMass = 0.0;
        
        // Find the best representative from each source in this cluster
        std::vector<std::pair<fluid::index, double>> sourceBests(nSources, {-1, 1e9});
        
        for (size_t idx : clusterMembers[c]) {
            const auto& m = allMasses[idx];
            double dist = std::abs(m.centerBin - clusterCenters[c]);
            
            if (dist < sourceBests[m.sourceIdx].second) {
                sourceBests[m.sourceIdx] = {m.massIdx, dist};
            }
            
            clusterMass += m.mass;
        }
        
        // Assign indices and check if all sources are represented
        for (size_t s = 0; s < nSources; ++s) {
            if (sourceBests[s].first != -1) {
                sourceIndices[s] = sourceBests[s].first;
            } else {
                // Find nearest mass from missing source
                double minDist = 1e9;
                fluid::index bestIdx = 0;
                
                for (size_t m = 0; m < masses[s].size(); ++m) {
                    double dist = std::abs(masses[s][m].centerBin - clusterCenters[c]);
                    if (dist < minDist) {
                        minDist = dist;
                        bestIdx = m;
                    }
                }
                
                sourceIndices[s] = bestIdx;
            }
        }
        
        // Add transport tuple
        transport.emplace_back(sourceIndices, clusterMass);
    }
    
    // Normalize transport masses
    double totalTransportMass = 0.0;
    for (const auto& t : transport) {
        totalTransportMass += std::get<1>(t);
    }
    
    if (totalTransportMass > 0) {
        for (auto& t : transport) {
            std::get<1>(t) /= totalTransportMass;
        }
    }
}

std::vector<AudioTransportN::SpetralMass> AudioTransportN::segmentSpectrum(
    const Eigen::Ref<Eigen::ArrayXd> mag,
    const Eigen::Ref<Eigen::ArrayXd> reassignedFreq,
    fluid::Allocator& alloc) {
    
    std::vector<SpetralMass> masses;
    double totalMass = mag.sum() + 1e-10;  // epsilon for stability
    
    // Compute sign changes in reassigned frequencies vs bin frequencies
    Eigen::ArrayXi sign = (reassignedFreq > mBinFreqs.head(reassignedFreq.size())).cast<int>();
    Eigen::ArrayXi changed = Eigen::ArrayXi::Zero(mBins);
    changed.segment(1, mBins - 1) = sign.segment(1, mBins - 1) - sign.segment(0, mBins - 1);
    
    SpetralMass currentMass{0, 0, 0, 0};
    
    for (fluid::index i = 1; i < mBins; i++) {
        if (changed(i) == -1) {
            // Local minimum found
            double d1 = reassignedFreq(i - 1) - mBinFreqs(i - 1);
            double d2 = mBinFreqs(i) - reassignedFreq(i);
            currentMass.centerBin = d1 < d2 ? i - 1 : i;
        }
        
        if (changed(i) == 1) {
            // End of current mass
            currentMass.endBin = i;
            currentMass.mass = mag.segment(currentMass.startBin, 
                                         i - currentMass.startBin).sum() / totalMass;
            masses.push_back(currentMass);
            currentMass = SpetralMass{i, i, i, 0};
        }
    }
    
    // Add final mass
    currentMass.endBin = mBins;
    currentMass.mass = mag.segment(currentMass.startBin, 
                                  mBins - currentMass.startBin).sum() / totalMass;
    masses.push_back(currentMass);
    
    return masses;
}

void AudioTransportN::placeMass(const SpetralMass mass, fluid::index bin, double scale,
                               double centerPhase, const Eigen::Ref<Eigen::ArrayXcd> input, 
                               Eigen::Ref<Eigen::ArrayXcd> output,
                               double nextPhase, Eigen::Ref<Eigen::ArrayXd> amplitudes, 
                               Eigen::Ref<Eigen::ArrayXd> phases) {
    
    double phaseShift = centerPhase - std::arg(input(mass.centerBin));
    
    for (fluid::index i = mass.startBin; i < mass.endBin; i++) {
        fluid::index pos = i + bin - mass.centerBin;
        
        if (pos < 0 || pos >= output.size()) continue;
        
        double phase = phaseShift + std::arg(input(i));
        double mag = scale * std::abs(input(i));
        
        output(pos) += std::polar(mag, phase);
        
        if (mag > amplitudes(pos)) {
            amplitudes(pos) = mag;
            phases(pos) = nextPhase;
        }
    }
}

} // namespace Explorer
} // namespace Acorex