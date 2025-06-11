#include "HDF5IO.h"
#include "Data.h"
#include "JSON.h"
#include <ofLog.h>
#include <filesystem>
#include <numeric>

namespace Utils {

// HDF5TransportIO implementation

std::string HDF5TransportIO::getFileGroupName(size_t fileIndex) const {
    return std::string(TRANSPORT_GROUP) + "/file_" + std::to_string(fileIndex);
}

std::string HDF5TransportIO::getFrameDatasetName(size_t fileIndex, const std::string& component) const {
    return getFileGroupName(fileIndex) + "/" + component;
}

void HDF5TransportIO::writeMetadata(H5::H5File& file, const TransportData& data, int compressionLevel) {
    try {
        // Create metadata group
        H5::Group metaGroup = file.createGroup(METADATA_GROUP);
        
        // Count total frames and get frame size
        size_t totalFrames = 0;
        size_t frameSize = 0;
        std::vector<size_t> frameCounts;
        
        for (size_t i = 0; i < data.fileCount(); ++i) {
            size_t count = data.frameCount(i);
            frameCounts.push_back(count);
            totalFrames += count;
            
            if (count > 0 && frameSize == 0) {
                const auto* frame = data.getFrame(i, 0);
                if (frame && frame->isValid()) {
                    frameSize = frame->frameSize();
                }
            }
        }
        
        // Write attributes
        H5::DataSpace scalar_space(H5S_SCALAR);
        
        H5::Attribute fileCountAttr = metaGroup.createAttribute(FILE_COUNT_ATTR, 
            H5::PredType::NATIVE_UINT64, scalar_space);
        size_t fileCount = data.fileCount();
        fileCountAttr.write(H5::PredType::NATIVE_UINT64, &fileCount);
        
        H5::Attribute frameCountAttr = metaGroup.createAttribute(FRAME_COUNT_ATTR, 
            H5::PredType::NATIVE_UINT64, scalar_space);
        frameCountAttr.write(H5::PredType::NATIVE_UINT64, &totalFrames);
        
        H5::Attribute frameSizeAttr = metaGroup.createAttribute(FRAME_SIZE_ATTR, 
            H5::PredType::NATIVE_UINT64, scalar_space);
        frameSizeAttr.write(H5::PredType::NATIVE_UINT64, &frameSize);
        
        H5::Attribute versionAttr = metaGroup.createAttribute(VERSION_ATTR, 
            H5::PredType::NATIVE_INT, scalar_space);
        int version = 1;
        versionAttr.write(H5::PredType::NATIVE_INT, &version);
        
        H5::Attribute compressionAttr = metaGroup.createAttribute(COMPRESSION_ATTR, 
            H5::PredType::NATIVE_INT, scalar_space);
        compressionAttr.write(H5::PredType::NATIVE_INT, &compressionLevel);
        
        // Write frame counts per file
        if (!frameCounts.empty()) {
            hsize_t dims[1] = {frameCounts.size()};
            H5::DataSpace frameCountsSpace(1, dims);
            H5::DataSet frameCountsDataset = metaGroup.createDataSet("frame_counts",
                H5::PredType::NATIVE_UINT64, frameCountsSpace);
            frameCountsDataset.write(frameCounts.data(), H5::PredType::NATIVE_UINT64);
        }
        
    } catch (H5::Exception& e) {
        ofLogError("HDF5TransportIO") << "Failed to write metadata: " << e.getCDetailMsg();
        throw;
    }
}

void HDF5TransportIO::readMetadata(H5::H5File& file) {
    try {
        H5::Group metaGroup = file.openGroup(METADATA_GROUP);
        
        // Read attributes
        H5::Attribute fileCountAttr = metaGroup.openAttribute(FILE_COUNT_ATTR);
        fileCountAttr.read(H5::PredType::NATIVE_UINT64, &m_fileCount);
        
        H5::Attribute frameSizeAttr = metaGroup.openAttribute(FRAME_SIZE_ATTR);
        frameSizeAttr.read(H5::PredType::NATIVE_UINT64, &m_frameSize);
        
        // Read frame counts
        if (metaGroup.nameExists("frame_counts")) {
            H5::DataSet frameCountsDataset = metaGroup.openDataSet("frame_counts");
            H5::DataSpace space = frameCountsDataset.getSpace();
            hsize_t dims[1];
            space.getSimpleExtentDims(dims);
            
            m_frameCounts.resize(dims[0]);
            frameCountsDataset.read(m_frameCounts.data(), H5::PredType::NATIVE_UINT64);
        }
        
    } catch (H5::Exception& e) {
        ofLogError("HDF5TransportIO") << "Failed to read metadata: " << e.getCDetailMsg();
        throw;
    }
}

H5::DataSet HDF5TransportIO::createCompressedDataset(H5::Group& group, const std::string& name,
                                                     const H5::DataSpace& space, int compressionLevel) {
    H5::DSetCreatPropList plist;
    
    if (compressionLevel > 0) {
        // Enable chunking (required for compression)
        hsize_t chunk_dims[2];
        int rank = space.getSimpleExtentNdims();
        hsize_t dims[2];
        space.getSimpleExtentDims(dims);
        
        // Set chunk size (balance between compression and access speed)
        if (rank == 2) {
            chunk_dims[0] = std::min(dims[0], hsize_t(100)); // Limit chunk to 100 frames
            chunk_dims[1] = dims[1]; // Full frame size
        } else {
            chunk_dims[0] = std::min(dims[0], hsize_t(100 * m_frameSize));
        }
        
        plist.setChunk(rank, chunk_dims);
        
        // Enable deflate (gzip) compression
        plist.setDeflate(compressionLevel);
        
        // Enable shuffle filter (improves compression for numeric data)
        plist.setShuffle();
    }
    
    return group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space, plist);
}

void HDF5TransportIO::writeFileData(H5::H5File& file, size_t fileIndex, 
                                   const std::vector<TransportFrame>& frames, int compressionLevel) {
    if (frames.empty()) {
        return;
    }
    
    try {
        // Create transport group if it doesn't exist
        H5::Group transportGroup;
        if (!file.nameExists(TRANSPORT_GROUP)) {
            transportGroup = file.createGroup(TRANSPORT_GROUP);
        } else {
            transportGroup = file.openGroup(TRANSPORT_GROUP);
        }
        
        // Create file group
        H5::Group fileGroup = transportGroup.createGroup("file_" + std::to_string(fileIndex));
        
        // Get dimensions
        size_t frameCount = frames.size();
        size_t frameSize = frames[0].frameSize();
        
        // Create 2D dataspace for each component
        hsize_t dims[2] = {frameCount, frameSize};
        H5::DataSpace dataSpace(2, dims);
        
        // Create datasets with compression
        H5::DataSet magnitudeDataset = createCompressedDataset(fileGroup, "magnitude", 
                                                              dataSpace, compressionLevel);
        H5::DataSet phaseDataset = createCompressedDataset(fileGroup, "phase", 
                                                           dataSpace, compressionLevel);
        H5::DataSet dHDataset = createCompressedDataset(fileGroup, "dH", 
                                                        dataSpace, compressionLevel);
        
        // Prepare contiguous data buffers
        std::vector<double> magnitudeData(frameCount * frameSize);
        std::vector<double> phaseData(frameCount * frameSize);
        std::vector<double> dHData(frameCount * frameSize);
        
        // Copy frame data to contiguous buffers
        for (size_t i = 0; i < frameCount; ++i) {
            const auto& frame = frames[i];
            if (frame.isValid() && frame.frameSize() == frameSize) {
                std::copy(frame.magnitude.begin(), frame.magnitude.end(),
                         magnitudeData.begin() + i * frameSize);
                std::copy(frame.phase.begin(), frame.phase.end(),
                         phaseData.begin() + i * frameSize);
                std::copy(frame.dH.begin(), frame.dH.end(),
                         dHData.begin() + i * frameSize);
            }
        }
        
        // Write data
        magnitudeDataset.write(magnitudeData.data(), H5::PredType::NATIVE_DOUBLE);
        phaseDataset.write(phaseData.data(), H5::PredType::NATIVE_DOUBLE);
        dHDataset.write(dHData.data(), H5::PredType::NATIVE_DOUBLE);
        
    } catch (H5::Exception& e) {
        ofLogError("HDF5TransportIO") << "Failed to write file data: " << e.getCDetailMsg();
        throw;
    }
}

void HDF5TransportIO::readFileData(H5::H5File& file, size_t fileIndex, 
                                  std::vector<TransportFrame>& frames) {
    try {
        std::string fileGroupPath = getFileGroupName(fileIndex);
        
        if (!file.nameExists(fileGroupPath)) {
            return;
        }
        
        H5::Group fileGroup = file.openGroup(fileGroupPath);
        
        // Open datasets
        H5::DataSet magnitudeDataset = fileGroup.openDataSet("magnitude");
        H5::DataSet phaseDataset = fileGroup.openDataSet("phase");
        H5::DataSet dHDataset = fileGroup.openDataSet("dH");
        
        // Get dimensions
        H5::DataSpace space = magnitudeDataset.getSpace();
        hsize_t dims[2];
        space.getSimpleExtentDims(dims);
        size_t frameCount = dims[0];
        size_t frameSize = dims[1];
        
        // Read data into contiguous buffers
        std::vector<double> magnitudeData(frameCount * frameSize);
        std::vector<double> phaseData(frameCount * frameSize);
        std::vector<double> dHData(frameCount * frameSize);
        
        magnitudeDataset.read(magnitudeData.data(), H5::PredType::NATIVE_DOUBLE);
        phaseDataset.read(phaseData.data(), H5::PredType::NATIVE_DOUBLE);
        dHDataset.read(dHData.data(), H5::PredType::NATIVE_DOUBLE);
        
        // Create frames
        frames.clear();
        frames.reserve(frameCount);
        
        for (size_t i = 0; i < frameCount; ++i) {
            TransportFrame frame(frameSize);
            
            std::copy(magnitudeData.begin() + i * frameSize,
                     magnitudeData.begin() + (i + 1) * frameSize,
                     frame.magnitude.begin());
            std::copy(phaseData.begin() + i * frameSize,
                     phaseData.begin() + (i + 1) * frameSize,
                     frame.phase.begin());
            std::copy(dHData.begin() + i * frameSize,
                     dHData.begin() + (i + 1) * frameSize,
                     frame.dH.begin());
            
            frames.push_back(std::move(frame));
        }
        
    } catch (H5::Exception& e) {
        ofLogError("HDF5TransportIO") << "Failed to read file data: " << e.getCDetailMsg();
        throw;
    }
}

bool HDF5TransportIO::write(const std::string& filepath, const TransportData& data, 
                           int compressionLevel) {
    try {
        // Disable HDF5 error printing to console
        H5::Exception::dontPrint();
        
        // Create HDF5 file
        m_file = std::make_unique<H5::H5File>(filepath, H5F_ACC_TRUNC);
        
        // Write metadata
        writeMetadata(*m_file, data, compressionLevel);
        
        // Track sizes for compression ratio
        size_t uncompressedSize = 0;
        
        // Write frame data for each file
        for (size_t i = 0; i < data.fileCount(); ++i) {
            if (data.frameCount(i) > 0) {
                writeFileData(*m_file, i, data.frames[i], compressionLevel);
                
                // Calculate uncompressed size
                for (const auto& frame : data.frames[i]) {
                    uncompressedSize += frame.frameSize() * sizeof(double) * 3;
                }
            }
        }
        
        // Close and reopen to get compressed size
        m_file->close();
        
        // Get compressed file size
        std::filesystem::path p(filepath);
        size_t compressedSize = std::filesystem::file_size(p);
        
        m_compressionRatio = uncompressedSize > 0 ? 
            static_cast<double>(uncompressedSize) / compressedSize : 1.0;
        
        ofLogNotice("HDF5TransportIO") << "Wrote HDF5 file with compression ratio: " 
                                       << m_compressionRatio << ":1";
        
        return true;
        
    } catch (H5::Exception& e) {
        ofLogError("HDF5TransportIO") << "HDF5 error: " << e.getCDetailMsg();
        return false;
    } catch (std::exception& e) {
        ofLogError("HDF5TransportIO") << "Error: " << e.what();
        return false;
    }
}

bool HDF5TransportIO::read(const std::string& filepath, TransportData& data) {
    try {
        if (!openForReading(filepath)) {
            return false;
        }
        
        data.clear();
        data.reserveFiles(m_fileCount);
        
        // Read all file data
        for (size_t i = 0; i < m_fileCount; ++i) {
            data.addFile();
            if (i < m_frameCounts.size() && m_frameCounts[i] > 0) {
                readFileData(*m_file, i, data.frames[i]);
            }
        }
        
        close();
        return true;
        
    } catch (H5::Exception& e) {
        ofLogError("HDF5TransportIO") << "HDF5 error: " << e.getCDetailMsg();
        return false;
    } catch (std::exception& e) {
        ofLogError("HDF5TransportIO") << "Error: " << e.what();
        return false;
    }
}

bool HDF5TransportIO::openForReading(const std::string& filepath) {
    try {
        close();
        
        // Disable HDF5 error printing to console
        H5::Exception::dontPrint();
        
        m_file = std::make_unique<H5::H5File>(filepath, H5F_ACC_RDONLY);
        
        // Read metadata
        readMetadata(*m_file);
        
        m_isOpen = true;
        return true;
        
    } catch (H5::Exception& e) {
        ofLogError("HDF5TransportIO") << "Failed to open HDF5 file: " << e.getCDetailMsg();
        return false;
    } catch (std::exception& e) {
        ofLogError("HDF5TransportIO") << "Error: " << e.what();
        return false;
    }
}

std::unique_ptr<TransportFrame> HDF5TransportIO::loadFrame(size_t fileIndex, size_t frameIndex) {
    if (!m_isOpen || !m_file) {
        return nullptr;
    }
    
    try {
        std::string fileGroupPath = getFileGroupName(fileIndex);
        
        if (!m_file->nameExists(fileGroupPath)) {
            return nullptr;
        }
        
        H5::Group fileGroup = m_file->openGroup(fileGroupPath);
        
        // Open datasets
        H5::DataSet magnitudeDataset = fileGroup.openDataSet("magnitude");
        H5::DataSet phaseDataset = fileGroup.openDataSet("phase");
        H5::DataSet dHDataset = fileGroup.openDataSet("dH");
        
        // Get dimensions
        H5::DataSpace fileSpace = magnitudeDataset.getSpace();
        hsize_t dims[2];
        fileSpace.getSimpleExtentDims(dims);
        
        if (frameIndex >= dims[0]) {
            return nullptr;
        }
        
        size_t frameSize = dims[1];
        
        // Define hyperslab for single frame
        hsize_t offset[2] = {frameIndex, 0};
        hsize_t count[2] = {1, frameSize};
        
        fileSpace.selectHyperslab(H5S_SELECT_SET, count, offset);
        
        // Define memory space
        hsize_t memDims[1] = {frameSize};
        H5::DataSpace memSpace(1, memDims);
        
        // Create frame
        auto frame = std::make_unique<TransportFrame>(frameSize);
        
        // Read data
        magnitudeDataset.read(frame->magnitude.data(), H5::PredType::NATIVE_DOUBLE, 
                             memSpace, fileSpace);
        
        phaseDataset.read(frame->phase.data(), H5::PredType::NATIVE_DOUBLE, 
                         memSpace, fileSpace);
        
        dHDataset.read(frame->dH.data(), H5::PredType::NATIVE_DOUBLE, 
                      memSpace, fileSpace);
        
        return frame;
        
    } catch (H5::Exception& e) {
        ofLogError("HDF5TransportIO") << "Failed to load frame: " << e.getCDetailMsg();
        return nullptr;
    }
}

bool HDF5TransportIO::getMetadata(size_t& fileCount, size_t& totalFrames) const {
    if (!m_isOpen) {
        return false;
    }
    
    fileCount = m_fileCount;
    totalFrames = std::accumulate(m_frameCounts.begin(), m_frameCounts.end(), size_t(0));
    return true;
}

double HDF5TransportIO::getCompressionRatio() const {
    return m_compressionRatio;
}

void HDF5TransportIO::close() {
    if (m_file) {
        try {
            m_file->close();
        } catch (H5::Exception& e) {
            // Ignore close errors
        }
        m_file.reset();
    }
    
    m_isOpen = false;
    m_fileCount = 0;
    m_frameSize = 0;
    m_frameCounts.clear();
    m_compressionRatio = 1.0;
}

// HDF5DataSetIO implementation

std::string HDF5DataSetIO::getHDF5Path(const std::string& jsonPath) {
    std::filesystem::path p(jsonPath);
    p.replace_extension(".h5");
    return p.string();
}

bool HDF5DataSetIO::write(const std::string& jsonPath, const DataSet& dataset,
                         int compressionLevel) {
    // Write main data to JSON (excluding transport data)
    DataSet datasetCopy = dataset;
    datasetCopy.transport.clear(); // Clear transport data from JSON
    
    if (!JSON::Write(jsonPath, datasetCopy)) {
        return false;
    }
    
    // Write transport data to HDF5 file if present
    if (dataset.analysisSettings.bTransport && dataset.transport.fileCount() > 0) {
        std::string hdf5Path = getHDF5Path(jsonPath);
        HDF5TransportIO io;
        
        if (!io.write(hdf5Path, dataset.transport, compressionLevel)) {
            // Clean up JSON file on failure
            std::filesystem::remove(jsonPath);
            return false;
        }
        
        ofLogNotice("HDF5DataSetIO") << "Compression ratio: " << io.getCompressionRatio() << ":1";
    }
    
    return true;
}

bool HDF5DataSetIO::read(const std::string& jsonPath, DataSet& dataset) {
    // Read main data from JSON
    if (!JSON::Read(jsonPath, dataset)) {
        return false;
    }
    
    // Read transport data from HDF5 file if expected
    if (dataset.analysisSettings.bTransport) {
        std::string hdf5Path = getHDF5Path(jsonPath);
        
        if (std::filesystem::exists(hdf5Path)) {
            HDF5TransportIO io;
            if (!io.read(hdf5Path, dataset.transport)) {
                // Continue without transport data, but log warning
                ofLogWarning("HDF5DataSetIO") << "Failed to load transport data from: " << hdf5Path;
                dataset.transport.clear();
            }
        }
    }
    
    return true;
}

} // namespace Utils