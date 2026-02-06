#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <cstring>
#include <stdexcept>
#include <cmath>
#include <limits>

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
    FLOAT16,  // new in v0.3.0 for ML/inference use cases
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
    if (dtype == "float16") return BandDataType::FLOAT16;
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
        case BandDataType::FLOAT16:
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

// Convert IEEE 754 half-precision (float16) to double
// Format: 1 sign bit, 5 exponent bits (bias 15), 10 mantissa bits
inline double float16_to_double(uint16_t h) {
    uint32_t sign = (h >> 15) & 0x1;
    uint32_t exp = (h >> 10) & 0x1F;
    uint32_t mant = h & 0x3FF;

    if (exp == 0) {
        if (mant == 0) {
            // Zero (positive or negative)
            return sign ? -0.0 : 0.0;
        } else {
            // Subnormal number
            double val = mant / 1024.0 * std::pow(2.0, -14);
            return sign ? -val : val;
        }
    } else if (exp == 31) {
        if (mant == 0) {
            // Infinity
            return sign ? -std::numeric_limits<double>::infinity()
                        : std::numeric_limits<double>::infinity();
        } else {
            // NaN
            return std::numeric_limits<double>::quiet_NaN();
        }
    } else {
        // Normalized number
        double val = (1.0 + mant / 1024.0) * std::pow(2.0, static_cast<int>(exp) - 15);
        return sign ? -val : val;
    }
}

// Decompress gzip data
std::vector<uint8_t> decompress_gzip(const uint8_t *data, size_t size);

// v0.4.0: Decompress JPEG image to raw RGB bytes
// Returns decoded image data in RGB format (3 bytes per pixel, row-major)
// width_out and height_out are set to the decoded image dimensions
std::vector<uint8_t> decompress_jpeg(const uint8_t *data, size_t size,
                                      int &width_out, int &height_out, int &channels_out);

// v0.4.0: Decompress WebP image to raw RGB/RGBA bytes
// Returns decoded image data in RGB or RGBA format (row-major)
// width_out and height_out are set to the decoded image dimensions
std::vector<uint8_t> decompress_webp(const uint8_t *data, size_t size,
                                      int &width_out, int &height_out, int &channels_out);

// Get single pixel value from raw (decompressed) band data
// data_size is the total size of the data buffer in bytes for bounds checking
inline double get_pixel_value(const uint8_t *data, size_t data_size, size_t offset, BandDataType dtype) {
    size_t elem_size = dtype_size(dtype);
    size_t byte_offset = offset * elem_size;
    if (byte_offset + elem_size > data_size) {
        throw std::out_of_range("Pixel read at byte offset " + std::to_string(byte_offset) +
                                " exceeds buffer of " + std::to_string(data_size) + " bytes");
    }
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
        case BandDataType::FLOAT16: {
            uint16_t val;
            std::memcpy(&val, data + offset * 2, 2);
            return float16_to_double(val);
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

// v0.4.0: Decode pixel from interleaved (BIP) layout
// For interleaved data, all bands are stored in a single 'pixels' column
// Data layout: [R0,G0,B0,R1,G1,B1,...] where each pixel's bands are adjacent
double decode_pixel_interleaved(const uint8_t *pixels_data, size_t pixels_size,
                                 const std::string &dtype_str,
                                 int pixel_x, int pixel_y, int width,
                                 int band_index, int num_bands,
                                 const std::string &compression);

// v0.4.0: Decode entire band from interleaved layout
std::vector<double> decode_band_interleaved(const uint8_t *pixels_data, size_t pixels_size,
                                             const std::string &dtype_str,
                                             int width, int height,
                                             int band_index, int num_bands,
                                             const std::string &compression);

// Statistics result structure
struct BandStats {
    int64_t count = 0;
    double sum = 0.0;
    double mean = 0.0;
    double min = 0.0;
    double max = 0.0;
    double stddev = 0.0;
};

// Compute statistics directly from band data without allocating full pixel array
// This is more memory efficient for large tiles
BandStats compute_band_stats(const uint8_t *band_data, size_t band_size,
                              const std::string &dtype_str,
                              int width, int height,
                              bool compressed,
                              bool has_nodata = false,
                              double nodata = 0.0);

} // namespace raquet
} // namespace duckdb
