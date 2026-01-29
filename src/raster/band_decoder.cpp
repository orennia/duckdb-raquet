#include "band_decoder.hpp"
#include <zlib.h>
#include <stdexcept>
#include <cstring>
#include <vector>
#include <string>
#include <limits>
#include <cmath>

namespace duckdb {
namespace raquet {

std::vector<uint8_t> decompress_gzip(const uint8_t *data, size_t size) {
    if (size == 0) {
        return {};
    }

    if (data == nullptr) {
        throw std::runtime_error("Null data pointer for decompression");
    }

    // For gzip: estimate decompressed size
    // A 256x256 tile with uint16 = 131072 bytes, float64 = 524288 bytes
    // Start with a larger initial estimate to avoid resizes
    size_t estimated_size = 256 * 256 * 8;  // Assume worst case (float64)
    if (size > 1000) {
        // For larger compressed data, use a ratio-based estimate
        estimated_size = std::max(estimated_size, size * 100);
    }

    std::vector<uint8_t> result(estimated_size);

    z_stream strm;
    memset(&strm, 0, sizeof(strm));
    strm.next_in = const_cast<Bytef*>(data);
    strm.avail_in = static_cast<uInt>(size);
    strm.next_out = result.data();
    strm.avail_out = static_cast<uInt>(result.size());

    // 15 + 32 = auto-detect gzip or zlib format
    int ret = inflateInit2(&strm, 15 + 32);
    if (ret != Z_OK) {
        throw std::runtime_error("inflateInit2 failed: " + std::to_string(ret));
    }

    // Use Z_NO_FLUSH and loop for potentially large data
    size_t total_out = 0;
    while (ret != Z_STREAM_END) {
        ret = inflate(&strm, Z_NO_FLUSH);

        if (ret == Z_BUF_ERROR || (ret == Z_OK && strm.avail_out == 0)) {
            // Need more output space
            total_out = result.size() - strm.avail_out;
            size_t new_size = result.size() * 2;
            result.resize(new_size);
            strm.next_out = result.data() + total_out;
            strm.avail_out = static_cast<uInt>(result.size() - total_out);
            continue;
        }

        if (ret != Z_OK && ret != Z_STREAM_END) {
            inflateEnd(&strm);
            throw std::runtime_error("inflate failed: " + std::to_string(ret));
        }
    }

    size_t decompressed_size = result.size() - strm.avail_out;
    inflateEnd(&strm);

    result.resize(decompressed_size);
    return result;
}

double decode_pixel(const uint8_t *band_data, size_t band_size,
                    const std::string &dtype_str,
                    int pixel_x, int pixel_y, int width,
                    bool compressed) {
    BandDataType dtype = parse_dtype(dtype_str);

    const uint8_t *data;
    std::vector<uint8_t> decompressed;

    if (compressed) {
        decompressed = decompress_gzip(band_data, band_size);
        data = decompressed.data();
    } else {
        data = band_data;
    }

    // Row-major order: offset = y * width + x
    size_t offset = static_cast<size_t>(pixel_y) * width + pixel_x;

    return get_pixel_value(data, offset, dtype);
}

std::vector<double> decode_band(const uint8_t *band_data, size_t band_size,
                                 const std::string &dtype_str,
                                 int width, int height,
                                 bool compressed) {
    BandDataType dtype = parse_dtype(dtype_str);

    const uint8_t *data;
    std::vector<uint8_t> decompressed;

    if (compressed) {
        decompressed = decompress_gzip(band_data, band_size);
        data = decompressed.data();
    } else {
        data = band_data;
    }

    size_t pixel_count = static_cast<size_t>(width) * height;
    std::vector<double> result(pixel_count);

    for (size_t i = 0; i < pixel_count; i++) {
        result[i] = get_pixel_value(data, i, dtype);
    }

    return result;
}

BandStats compute_band_stats(const uint8_t *band_data, size_t band_size,
                              const std::string &dtype_str,
                              int width, int height,
                              bool compressed,
                              bool has_nodata,
                              double nodata) {
    BandDataType dtype = parse_dtype(dtype_str);

    const uint8_t *data;
    std::vector<uint8_t> decompressed;

    if (compressed) {
        decompressed = decompress_gzip(band_data, band_size);
        data = decompressed.data();
    } else {
        data = band_data;
    }

    size_t pixel_count = static_cast<size_t>(width) * height;

    BandStats stats;
    stats.min = std::numeric_limits<double>::max();
    stats.max = std::numeric_limits<double>::lowest();

    // Welford's online algorithm variables
    double m2 = 0.0;

    for (size_t i = 0; i < pixel_count; i++) {
        double val = get_pixel_value(data, i, dtype);

        // Skip nodata values (handle NaN specially since NaN != NaN)
        if (has_nodata && (val == nodata || (std::isnan(val) && std::isnan(nodata)))) {
            continue;
        }

        stats.count++;
        stats.sum += val;

        if (val < stats.min) stats.min = val;
        if (val > stats.max) stats.max = val;

        // Welford's algorithm for mean and variance
        double delta = val - stats.mean;
        stats.mean += delta / stats.count;
        double delta2 = val - stats.mean;
        m2 += delta * delta2;
    }

    // Finalize statistics
    if (stats.count > 0) {
        stats.mean = stats.sum / stats.count;
        if (stats.count > 1) {
            stats.stddev = std::sqrt(m2 / (stats.count - 1));
        }
    } else {
        stats.min = 0;
        stats.max = 0;
    }

    return stats;
}

} // namespace raquet
} // namespace duckdb
