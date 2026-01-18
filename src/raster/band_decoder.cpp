#include "band_decoder.hpp"
#include <zlib.h>
#include <stdexcept>

namespace duckdb {
namespace raquet {

std::vector<uint8_t> decompress_gzip(const uint8_t *data, size_t size) {
    if (size == 0) {
        return {};
    }

    // Initialize zlib for gzip decompression
    z_stream stream = {};
    stream.zalloc = Z_NULL;
    stream.zfree = Z_NULL;
    stream.opaque = Z_NULL;
    stream.avail_in = static_cast<uInt>(size);
    stream.next_in = const_cast<Bytef*>(data);

    // 15 + 16 for gzip format, 15 + 32 for auto-detect gzip/zlib
    int init_ret = inflateInit2(&stream, 15 + 32);
    if (init_ret != Z_OK) {
        throw std::runtime_error("Failed to initialize zlib for decompression");
    }

    std::vector<uint8_t> result;
    result.reserve(size * 4);  // Initial estimate

    uint8_t buffer[32768];
    int ret;
    size_t max_iterations = (size * 100) / sizeof(buffer) + 100;  // Safety limit
    size_t iterations = 0;

    do {
        stream.avail_out = sizeof(buffer);
        stream.next_out = buffer;

        ret = inflate(&stream, Z_NO_FLUSH);

        if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR || ret == Z_NEED_DICT) {
            inflateEnd(&stream);
            throw std::runtime_error("Zlib decompression error: " + std::to_string(ret));
        }

        size_t have = sizeof(buffer) - stream.avail_out;
        if (have > 0) {
            result.insert(result.end(), buffer, buffer + have);
        }

        iterations++;
        if (iterations > max_iterations) {
            inflateEnd(&stream);
            throw std::runtime_error("Decompression exceeded maximum iterations");
        }

    } while (ret != Z_STREAM_END && stream.avail_in > 0);

    inflateEnd(&stream);
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

} // namespace raquet
} // namespace duckdb
