#pragma once

#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <H5Cpp.h>
#include "Data.h"

namespace Acorex {
namespace Utils {

// Forward declarations are no longer needed since we include Data.h

// HDF5-based TransportData I/O
class HDF5TransportIO {
public:
    HDF5TransportIO() = default;
    ~HDF5TransportIO() = default;
    
    // Write TransportData to HDF5 file with compression
    bool write(const std::string& filepath, const TransportData& data, 
               int compressionLevel = 6); // 0-9, where 0 is no compression
    
    // Read TransportData from HDF5 file
    bool read(const std::string& filepath, TransportData& data);
    
    // Lazy loading support - open file without loading all data
    bool openForReading(const std::string& filepath);
    
    // Load specific frame on demand
    std::unique_ptr<TransportFrame> loadFrame(size_t fileIndex, size_t frameIndex);
    
    // Get metadata without loading data
    bool getMetadata(size_t& fileCount, size_t& totalFrames) const;
    
    // Get compression ratio info
    double getCompressionRatio() const;
    
    // Close the file
    void close();
    
private:
    // HDF5 dataset names
    static constexpr const char* TRANSPORT_GROUP = "/transport";
    static constexpr const char* METADATA_GROUP = "/metadata";
    static constexpr const char* FILE_COUNT_ATTR = "file_count";
    static constexpr const char* FRAME_COUNT_ATTR = "frame_count";
    static constexpr const char* FRAME_SIZE_ATTR = "frame_size";
    static constexpr const char* VERSION_ATTR = "version";
    static constexpr const char* COMPRESSION_ATTR = "compression_level";
    
    // Helper methods
    std::string getFileGroupName(size_t fileIndex) const;
    std::string getFrameDatasetName(size_t fileIndex, const std::string& component) const;
    
    void writeMetadata(H5::H5File& file, const TransportData& data, int compressionLevel);
    void readMetadata(H5::H5File& file);
    
    void writeFileData(H5::H5File& file, size_t fileIndex, 
                      const std::vector<TransportFrame>& frames, int compressionLevel);
    void readFileData(H5::H5File& file, size_t fileIndex, 
                     std::vector<TransportFrame>& frames);
    
    // Create dataset with compression
    H5::DataSet createCompressedDataset(H5::Group& group, const std::string& name,
                                       const H5::DataSpace& space, int compressionLevel);
    
    // Member variables
    std::unique_ptr<H5::H5File> m_file;
    size_t m_fileCount = 0;
    size_t m_frameSize = 0;
    std::vector<size_t> m_frameCounts;
    bool m_isOpen = false;
    double m_compressionRatio = 1.0;
};

// HDF5-based DataSet I/O with hybrid approach
class HDF5DataSetIO {
public:
    // Write DataSet with transport data in HDF5 format
    static bool write(const std::string& jsonPath, const DataSet& dataset,
                     int compressionLevel = 6);
    
    // Read DataSet with HDF5 transport data
    static bool read(const std::string& jsonPath, DataSet& dataset);
    
private:
    // Generate HDF5 filename from JSON filename
    static std::string getHDF5Path(const std::string& jsonPath);
};

} // namespace Utils
} // namespace Acorex