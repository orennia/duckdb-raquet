#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <cstring>
#include <stdexcept>

namespace duckdb {
namespace raquet {

// Supported data types for band pixels
enum class BandDataType {
    UINT8,
    INT8,
    UINT16,
    INT16,
    UINT32,
    INT32,
    UINT64,
    INT64,
    FLOAT32,
    FLOAT64
};

// Parse data type from string
inline BandDataType parse_dtype(const std::string &dtype) {
    if (dtype == "uint8") return BandDataType::UINT8;
    if (dtype == "int8") return BandDataType::INT8;
    if (dtype == "uint16") return BandDataType::UINT16;
    if (dtype == "int16") return BandDataType::INT16;
    if (dtype == "uint32") return BandDataType::UINT32;
    if (dtype == "int32") return BandDataType::INT32;
    if (dtype == "uint64") return BandDataType::UINT64;
    if (dtype == "int64") return BandDataType::INT64;
    if (dtype == "float32") return BandDataType::FLOAT32;
    if (dtype == "float64") return BandDataType::FLOAT64;
    throw std::invalid_argument("Unknown data type: " + dtype);
}

// Get byte size for data type
inline size_t dtype_size(BandDataType dtype) {
    switch (dtype) {
        case BandDataType::UINT8:
        case BandDataType::INT8:
            return 1;
        case BandDataType::UINT16:
        case BandDataType::INT16:
            return 2;
        case BandDataType::UINT32:
        case BandDataType::INT32:
        case BandDataType::FLOAT32:
            return 4;
        case BandDataType::UINT64:
        case BandDataType::INT64:
        case BandDataType::FLOAT64:
            return 8;
    }
    return 0;
}

// Decompress gzip data
std::vector<uint8_t> decompress_gzip(const uint8_t *data, size_t size);

// Get single pixel value from raw (decompressed) band data
inline double get_pixel_value(const uint8_t *data, size_t offset, BandDataType dtype) {
    switch (dtype) {
        case BandDataType::UINT8:
            return static_cast<double>(data[offset]);
        case BandDataType::INT8:
            return static_cast<double>(reinterpret_cast<const int8_t*>(data)[offset]);
        case BandDataType::UINT16: {
            uint16_t val;
            std::memcpy(&val, data + offset * 2, 2);
            return static_cast<double>(val);
        }
        case BandDataType::INT16: {
            int16_t val;
            std::memcpy(&val, data + offset * 2, 2);
            return static_cast<double>(val);
        }
        case BandDataType::UINT32: {
            uint32_t val;
            std::memcpy(&val, data + offset * 4, 4);
            return static_cast<double>(val);
        }
        case BandDataType::INT32: {
            int32_t val;
            std::memcpy(&val, data + offset * 4, 4);
            return static_cast<double>(val);
        }
        case BandDataType::UINT64: {
            uint64_t val;
            std::memcpy(&val, data + offset * 8, 8);
            return static_cast<double>(val);
        }
        case BandDataType::INT64: {
            int64_t val;
            std::memcpy(&val, data + offset * 8, 8);
            return static_cast<double>(val);
        }
        case BandDataType::FLOAT32: {
            float val;
            std::memcpy(&val, data + offset * 4, 4);
            return static_cast<double>(val);
        }
        case BandDataType::FLOAT64: {
            double val;
            std::memcpy(&val, data + offset * 8, 8);
            return val;
        }
    }
    return 0.0;
}

// Decode band data and get pixel value at x, y
double decode_pixel(const uint8_t *band_data, size_t band_size,
                    const std::string &dtype_str,
                    int pixel_x, int pixel_y, int width,
                    bool compressed);

// Decode entire band to vector of doubles
std::vector<double> decode_band(const uint8_t *band_data, size_t band_size,
                                 const std::string &dtype_str,
                                 int width, int height,
                                 bool compressed);

} // namespace raquet
} // namespace duckdb
