#pragma once

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <stdexcept>
#include <filesystem>

#ifdef _WIN32
    #include <windows.h>
#else
    #include <sys/mman.h>
    #include <sys/stat.h>
    #include <fcntl.h>
    #include <unistd.h>
#endif

namespace Utils {

// Forward declarations
struct TransportFrame;
struct TransportData;

// Memory-mapped file wrapper class
class MemoryMappedFile {
public:
    MemoryMappedFile() = default;
    ~MemoryMappedFile();
    
    // Disable copy operations
    MemoryMappedFile(const MemoryMappedFile&) = delete;
    MemoryMappedFile& operator=(const MemoryMappedFile&) = delete;
    
    // Enable move operations
    MemoryMappedFile(MemoryMappedFile&& other) noexcept;
    MemoryMappedFile& operator=(MemoryMappedFile&& other) noexcept;
    
    // Open file for reading or writing
    bool open(const std::string& filepath, bool readOnly = true);
    
    // Create new file with specified size
    bool create(const std::string& filepath, size_t size);
    
    // Close the mapped file
    void close();
    
    // Get raw data pointer
    void* data() { return m_data; }
    const void* data() const { return m_data; }
    
    // Get file size
    size_t size() const { return m_size; }
    
    // Check if file is open
    bool isOpen() const { return m_data != nullptr; }
    
private:
    void* m_data = nullptr;
    size_t m_size = 0;
    
#ifdef _WIN32
    HANDLE m_fileHandle = INVALID_HANDLE_VALUE;
    HANDLE m_mapHandle = nullptr;
#else
    int m_fd = -1;
#endif
};

// Binary format structure for TransportFrame storage
struct TransportFrameHeader {
    uint32_t magic = 0x54524654; // 'TRFT'
    uint32_t version = 1;
    uint32_t frameSize = 0;
    uint32_t flags = 0;
};

// Binary format structure for TransportData file
struct TransportDataHeader {
    uint32_t magic = 0x54524454; // 'TRDT'
    uint32_t version = 1;
    uint32_t fileCount = 0;
    uint32_t totalFrames = 0;
    uint64_t indexTableOffset = 0;
    uint64_t dataOffset = 0;
};

// Index entry for fast frame lookup
struct FrameIndexEntry {
    uint32_t fileIndex;
    uint32_t frameIndex;
    uint64_t offset;
    uint32_t size;
};

// Memory-mapped TransportData I/O
class TransportDataIO {
public:
    TransportDataIO() = default;
    ~TransportDataIO() = default;
    
    // Write TransportData to memory-mapped file
    bool write(const std::string& filepath, const TransportData& data);
    
    // Read TransportData from memory-mapped file
    bool read(const std::string& filepath, TransportData& data);
    
    // Lazy loading support - open file without loading all data
    bool openForReading(const std::string& filepath);
    
    // Load specific frame on demand
    std::unique_ptr<TransportFrame> loadFrame(size_t fileIndex, size_t frameIndex);
    
    // Get metadata without loading data
    bool getMetadata(uint32_t& fileCount, uint32_t& totalFrames) const;
    
    // Close the file
    void close();
    
private:
    // Helper methods for serialization
    size_t calculateFileSize(const TransportData& data) const;
    size_t serializeFrame(void* dest, const TransportFrame& frame) const;
    bool deserializeFrame(const void* src, size_t size, TransportFrame& frame) const;
    
    // Build index table for fast lookup
    void buildIndexTable(const TransportData& data, std::vector<FrameIndexEntry>& index) const;
    
    // Member variables
    MemoryMappedFile m_file;
    TransportDataHeader m_header;
    std::vector<FrameIndexEntry> m_indexTable;
    bool m_isOpen = false;
};

// Hybrid approach: JSON metadata + memory-mapped transport data
class HybridDataSetIO {
public:
    // Write DataSet with transport data in separate memory-mapped file
    static bool write(const std::string& jsonPath, const DataSet& dataset);
    
    // Read DataSet with lazy-loaded transport data
    static bool read(const std::string& jsonPath, DataSet& dataset);
    
private:
    // Generate binary filename from JSON filename
    static std::string getBinaryPath(const std::string& jsonPath);
};

} // namespace Utils