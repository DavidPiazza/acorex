#include "MemoryMappedIO.h"
#include "Data.h"
#include "JSON.h"
#include <ofLog.h>
#include <cstring>
#include <algorithm>
#include <fstream>

namespace Acorex {
namespace Utils {

// MemoryMappedFile implementation

MemoryMappedFile::~MemoryMappedFile() {
    close();
}

MemoryMappedFile::MemoryMappedFile(MemoryMappedFile&& other) noexcept
    : m_data(other.m_data)
    , m_size(other.m_size)
#ifdef _WIN32
    , m_fileHandle(other.m_fileHandle)
    , m_mapHandle(other.m_mapHandle)
#else
    , m_fd(other.m_fd)
#endif
{
    other.m_data = nullptr;
    other.m_size = 0;
#ifdef _WIN32
    other.m_fileHandle = INVALID_HANDLE_VALUE;
    other.m_mapHandle = nullptr;
#else
    other.m_fd = -1;
#endif
}

MemoryMappedFile& MemoryMappedFile::operator=(MemoryMappedFile&& other) noexcept {
    if (this != &other) {
        close();
        
        m_data = other.m_data;
        m_size = other.m_size;
#ifdef _WIN32
        m_fileHandle = other.m_fileHandle;
        m_mapHandle = other.m_mapHandle;
#else
        m_fd = other.m_fd;
#endif
        
        other.m_data = nullptr;
        other.m_size = 0;
#ifdef _WIN32
        other.m_fileHandle = INVALID_HANDLE_VALUE;
        other.m_mapHandle = nullptr;
#else
        other.m_fd = -1;
#endif
    }
    return *this;
}

bool MemoryMappedFile::open(const std::string& filepath, bool readOnly) {
    close();
    
#ifdef _WIN32
    DWORD access = readOnly ? GENERIC_READ : (GENERIC_READ | GENERIC_WRITE);
    DWORD shareMode = readOnly ? FILE_SHARE_READ : 0;
    DWORD createMode = readOnly ? OPEN_EXISTING : OPEN_ALWAYS;
    
    m_fileHandle = CreateFileA(filepath.c_str(), access, shareMode, nullptr, createMode, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (m_fileHandle == INVALID_HANDLE_VALUE) {
        return false;
    }
    
    LARGE_INTEGER fileSize;
    if (!GetFileSizeEx(m_fileHandle, &fileSize)) {
        CloseHandle(m_fileHandle);
        m_fileHandle = INVALID_HANDLE_VALUE;
        return false;
    }
    
    m_size = static_cast<size_t>(fileSize.QuadPart);
    
    if (m_size == 0) {
        return true; // Empty file, but successfully opened
    }
    
    DWORD protect = readOnly ? PAGE_READONLY : PAGE_READWRITE;
    m_mapHandle = CreateFileMappingA(m_fileHandle, nullptr, protect, 0, 0, nullptr);
    if (!m_mapHandle) {
        CloseHandle(m_fileHandle);
        m_fileHandle = INVALID_HANDLE_VALUE;
        return false;
    }
    
    DWORD mapAccess = readOnly ? FILE_MAP_READ : FILE_MAP_ALL_ACCESS;
    m_data = MapViewOfFile(m_mapHandle, mapAccess, 0, 0, 0);
    if (!m_data) {
        CloseHandle(m_mapHandle);
        CloseHandle(m_fileHandle);
        m_mapHandle = nullptr;
        m_fileHandle = INVALID_HANDLE_VALUE;
        return false;
    }
#else
    int flags = readOnly ? O_RDONLY : O_RDWR;
    m_fd = ::open(filepath.c_str(), flags);
    if (m_fd < 0) {
        return false;
    }
    
    struct stat st;
    if (fstat(m_fd, &st) != 0) {
        ::close(m_fd);
        m_fd = -1;
        return false;
    }
    
    m_size = static_cast<size_t>(st.st_size);
    
    if (m_size == 0) {
        return true; // Empty file, but successfully opened
    }
    
    int prot = readOnly ? PROT_READ : (PROT_READ | PROT_WRITE);
    m_data = mmap(nullptr, m_size, prot, MAP_SHARED, m_fd, 0);
    if (m_data == MAP_FAILED) {
        ::close(m_fd);
        m_fd = -1;
        m_data = nullptr;
        return false;
    }
#endif
    
    return true;
}

bool MemoryMappedFile::create(const std::string& filepath, size_t size) {
    close();
    
    if (size == 0) {
        return false;
    }
    
#ifdef _WIN32
    m_fileHandle = CreateFileA(filepath.c_str(), GENERIC_READ | GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (m_fileHandle == INVALID_HANDLE_VALUE) {
        return false;
    }
    
    LARGE_INTEGER fileSize;
    fileSize.QuadPart = static_cast<LONGLONG>(size);
    if (!SetFilePointerEx(m_fileHandle, fileSize, nullptr, FILE_BEGIN) || !SetEndOfFile(m_fileHandle)) {
        CloseHandle(m_fileHandle);
        m_fileHandle = INVALID_HANDLE_VALUE;
        return false;
    }
    
    m_size = size;
    
    m_mapHandle = CreateFileMappingA(m_fileHandle, nullptr, PAGE_READWRITE, 0, 0, nullptr);
    if (!m_mapHandle) {
        CloseHandle(m_fileHandle);
        m_fileHandle = INVALID_HANDLE_VALUE;
        return false;
    }
    
    m_data = MapViewOfFile(m_mapHandle, FILE_MAP_ALL_ACCESS, 0, 0, 0);
    if (!m_data) {
        CloseHandle(m_mapHandle);
        CloseHandle(m_fileHandle);
        m_mapHandle = nullptr;
        m_fileHandle = INVALID_HANDLE_VALUE;
        return false;
    }
#else
    m_fd = ::open(filepath.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (m_fd < 0) {
        return false;
    }
    
    if (ftruncate(m_fd, static_cast<off_t>(size)) != 0) {
        ::close(m_fd);
        m_fd = -1;
        return false;
    }
    
    m_size = size;
    
    m_data = mmap(nullptr, m_size, PROT_READ | PROT_WRITE, MAP_SHARED, m_fd, 0);
    if (m_data == MAP_FAILED) {
        ::close(m_fd);
        m_fd = -1;
        m_data = nullptr;
        return false;
    }
#endif
    
    return true;
}

void MemoryMappedFile::close() {
    if (m_data) {
#ifdef _WIN32
        UnmapViewOfFile(m_data);
#else
        munmap(m_data, m_size);
#endif
        m_data = nullptr;
    }
    
#ifdef _WIN32
    if (m_mapHandle) {
        CloseHandle(m_mapHandle);
        m_mapHandle = nullptr;
    }
    
    if (m_fileHandle != INVALID_HANDLE_VALUE) {
        CloseHandle(m_fileHandle);
        m_fileHandle = INVALID_HANDLE_VALUE;
    }
#else
    if (m_fd >= 0) {
        ::close(m_fd);
        m_fd = -1;
    }
#endif
    
    m_size = 0;
}

// TransportDataIO implementation

size_t TransportDataIO::calculateFileSize(const TransportData& data) const {
    size_t totalSize = sizeof(TransportDataHeader);
    
    // Calculate index table size
    size_t indexEntries = 0;
    for (size_t i = 0; i < data.fileCount(); ++i) {
        indexEntries += data.frameCount(i);
    }
    totalSize += indexEntries * sizeof(FrameIndexEntry);
    
    // Calculate data size
    for (size_t i = 0; i < data.fileCount(); ++i) {
        for (size_t j = 0; j < data.frameCount(i); ++j) {
            const auto* frame = data.getFrame(i, j);
            if (frame && frame->isValid()) {
                // Each frame: header + 3 vectors of doubles
                totalSize += sizeof(TransportFrameHeader);
                totalSize += frame->frameSize() * sizeof(double) * 3;
            }
        }
    }
    
    return totalSize;
}

size_t TransportDataIO::serializeFrame(void* dest, const TransportFrame& frame) const {
    if (!frame.isValid()) {
        return 0;
    }
    
    uint8_t* ptr = static_cast<uint8_t*>(dest);
    
    // Write header
    TransportFrameHeader header;
    header.frameSize = static_cast<uint32_t>(frame.frameSize());
    std::memcpy(ptr, &header, sizeof(header));
    ptr += sizeof(header);
    
    // Write magnitude
    size_t vectorSize = frame.frameSize() * sizeof(double);
    std::memcpy(ptr, frame.magnitude.data(), vectorSize);
    ptr += vectorSize;
    
    // Write phase
    std::memcpy(ptr, frame.phase.data(), vectorSize);
    ptr += vectorSize;
    
    // Write dH
    std::memcpy(ptr, frame.dH.data(), vectorSize);
    ptr += vectorSize;
    
    return sizeof(header) + 3 * vectorSize;
}

bool TransportDataIO::deserializeFrame(const void* src, size_t size, TransportFrame& frame) const {
    if (size < sizeof(TransportFrameHeader)) {
        return false;
    }
    
    const uint8_t* ptr = static_cast<const uint8_t*>(src);
    
    // Read header
    TransportFrameHeader header;
    std::memcpy(&header, ptr, sizeof(header));
    ptr += sizeof(header);
    
    if (header.magic != 0x54524654 || header.version != 1) {
        return false;
    }
    
    size_t expectedSize = sizeof(header) + header.frameSize * sizeof(double) * 3;
    if (size < expectedSize) {
        return false;
    }
    
    // Read vectors
    frame.resize(header.frameSize);
    
    size_t vectorSize = header.frameSize * sizeof(double);
    std::memcpy(frame.magnitude.data(), ptr, vectorSize);
    ptr += vectorSize;
    
    std::memcpy(frame.phase.data(), ptr, vectorSize);
    ptr += vectorSize;
    
    std::memcpy(frame.dH.data(), ptr, vectorSize);
    
    return true;
}

void TransportDataIO::buildIndexTable(const TransportData& data, std::vector<FrameIndexEntry>& index) const {
    index.clear();
    
    uint64_t currentOffset = sizeof(TransportDataHeader);
    
    // Reserve space for index table
    size_t totalFrames = 0;
    for (size_t i = 0; i < data.fileCount(); ++i) {
        totalFrames += data.frameCount(i);
    }
    
    currentOffset += totalFrames * sizeof(FrameIndexEntry);
    
    // Build index entries
    for (size_t i = 0; i < data.fileCount(); ++i) {
        for (size_t j = 0; j < data.frameCount(i); ++j) {
            const auto* frame = data.getFrame(i, j);
            if (frame && frame->isValid()) {
                FrameIndexEntry entry;
                entry.fileIndex = static_cast<uint32_t>(i);
                entry.frameIndex = static_cast<uint32_t>(j);
                entry.offset = currentOffset;
                
                // Calculate frame size
                size_t frameSize = sizeof(TransportFrameHeader) + frame->frameSize() * sizeof(double) * 3;
                entry.size = static_cast<uint32_t>(frameSize);
                
                index.push_back(entry);
                currentOffset += frameSize;
            }
        }
    }
}

bool TransportDataIO::write(const std::string& filepath, const TransportData& data) {
    // Calculate required file size
    size_t fileSize = calculateFileSize(data);
    
    // Create memory-mapped file
    if (!m_file.create(filepath, fileSize)) {
        return false;
    }
    
    uint8_t* ptr = static_cast<uint8_t*>(m_file.data());
    
    // Build index table
    std::vector<FrameIndexEntry> indexTable;
    buildIndexTable(data, indexTable);
    
    // Write header
    TransportDataHeader header;
    header.fileCount = static_cast<uint32_t>(data.fileCount());
    header.totalFrames = static_cast<uint32_t>(indexTable.size());
    header.indexTableOffset = sizeof(TransportDataHeader);
    header.dataOffset = header.indexTableOffset + indexTable.size() * sizeof(FrameIndexEntry);
    
    std::memcpy(ptr, &header, sizeof(header));
    ptr += sizeof(header);
    
    // Write index table
    if (!indexTable.empty()) {
        std::memcpy(ptr, indexTable.data(), indexTable.size() * sizeof(FrameIndexEntry));
        ptr += indexTable.size() * sizeof(FrameIndexEntry);
    }
    
    // Write frame data
    for (const auto& entry : indexTable) {
        const auto* frame = data.getFrame(entry.fileIndex, entry.frameIndex);
        if (frame) {
            size_t written = serializeFrame(ptr, *frame);
            ptr += written;
        }
    }
    
    m_file.close();
    return true;
}

bool TransportDataIO::read(const std::string& filepath, TransportData& data) {
    if (!openForReading(filepath)) {
        return false;
    }
    
    data.clear();
    data.reserveFiles(m_header.fileCount);
    
    // Create file structure
    std::vector<uint32_t> frameCounts(m_header.fileCount, 0);
    for (const auto& entry : m_indexTable) {
        frameCounts[entry.fileIndex] = std::max(frameCounts[entry.fileIndex], entry.frameIndex + 1);
    }
    
    for (uint32_t i = 0; i < m_header.fileCount; ++i) {
        data.addFile();
        data.reserveFrames(i, frameCounts[i]);
        
        // Initialize with empty frames
        for (uint32_t j = 0; j < frameCounts[i]; ++j) {
            data.frames[i].emplace_back();
        }
    }
    
    // Load all frames
    const uint8_t* basePtr = static_cast<const uint8_t*>(m_file.data());
    for (const auto& entry : m_indexTable) {
        TransportFrame frame;
        if (deserializeFrame(basePtr + entry.offset, entry.size, frame)) {
            data.frames[entry.fileIndex][entry.frameIndex] = std::move(frame);
        }
    }
    
    close();
    return true;
}

bool TransportDataIO::openForReading(const std::string& filepath) {
    close();
    
    if (!m_file.open(filepath, true)) {
        return false;
    }
    
    if (m_file.size() < sizeof(TransportDataHeader)) {
        m_file.close();
        return false;
    }
    
    // Read header
    std::memcpy(&m_header, m_file.data(), sizeof(m_header));
    
    if (m_header.magic != 0x54524454 || m_header.version != 1) {
        m_file.close();
        return false;
    }
    
    // Read index table
    if (m_header.totalFrames > 0) {
        m_indexTable.resize(m_header.totalFrames);
        const uint8_t* indexPtr = static_cast<const uint8_t*>(m_file.data()) + m_header.indexTableOffset;
        std::memcpy(m_indexTable.data(), indexPtr, m_header.totalFrames * sizeof(FrameIndexEntry));
    }
    
    m_isOpen = true;
    return true;
}

std::unique_ptr<TransportFrame> TransportDataIO::loadFrame(size_t fileIndex, size_t frameIndex) {
    if (!m_isOpen) {
        return nullptr;
    }
    
    // Find frame in index
    auto it = std::find_if(m_indexTable.begin(), m_indexTable.end(),
        [fileIndex, frameIndex](const FrameIndexEntry& entry) {
            return entry.fileIndex == fileIndex && entry.frameIndex == frameIndex;
        });
    
    if (it == m_indexTable.end()) {
        return nullptr;
    }
    
    auto frame = std::make_unique<TransportFrame>();
    const uint8_t* ptr = static_cast<const uint8_t*>(m_file.data()) + it->offset;
    
    if (deserializeFrame(ptr, it->size, *frame)) {
        return frame;
    }
    
    return nullptr;
}

bool TransportDataIO::getMetadata(uint32_t& fileCount, uint32_t& totalFrames) const {
    if (!m_isOpen) {
        return false;
    }
    
    fileCount = m_header.fileCount;
    totalFrames = m_header.totalFrames;
    return true;
}

void TransportDataIO::close() {
    m_file.close();
    m_indexTable.clear();
    m_isOpen = false;
    m_header = TransportDataHeader();
}

// HybridDataSetIO implementation

std::string HybridDataSetIO::getBinaryPath(const std::string& jsonPath) {
    std::filesystem::path p(jsonPath);
    p.replace_extension(".transport");
    return p.string();
}

bool HybridDataSetIO::write(const std::string& jsonPath, const DataSet& dataset) {
    // Write main data to JSON (excluding transport data)
    DataSet datasetCopy = dataset;
    datasetCopy.transport.clear(); // Clear transport data from JSON
    
    JSON json;
    if (!json.Write(jsonPath, datasetCopy)) {
        return false;
    }
    
    // Write transport data to binary file if present
    if (dataset.analysisSettings.bTransport && dataset.transport.fileCount() > 0) {
        std::string binaryPath = getBinaryPath(jsonPath);
        TransportDataIO io;
        
        if (!io.write(binaryPath, dataset.transport)) {
            // Clean up JSON file on failure
            std::filesystem::remove(jsonPath);
            return false;
        }
    }
    
    return true;
}

bool HybridDataSetIO::read(const std::string& jsonPath, DataSet& dataset) {
    // Read main data from JSON
    JSON json;
    if (!json.Read(jsonPath, dataset)) {
        return false;
    }
    
    // Read transport data from binary file if expected
    if (dataset.analysisSettings.bTransport) {
        std::string binaryPath = getBinaryPath(jsonPath);
        
        if (std::filesystem::exists(binaryPath)) {
            TransportDataIO io;
            if (!io.read(binaryPath, dataset.transport)) {
                // Continue without transport data, but log warning
                ofLogWarning("HybridDataSetIO") << "Failed to load transport data from: " << binaryPath;
                dataset.transport.clear();
            }
        }
    }
    
    return true;
}

} // namespace Utils
} // namespace Acorex