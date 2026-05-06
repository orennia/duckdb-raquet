#pragma once

#include "band_decoder.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {
namespace raquet {

// Compress raw data with gzip
std::vector<uint8_t> compress_gzip(const uint8_t *data, size_t size);

// Encode raw RGB/grayscale pixels as JPEG
// Input: row-major pixel data, width * height * channels bytes
// Returns JPEG-encoded bytes
std::vector<uint8_t> encode_jpeg(const uint8_t *data, int width, int height,
                                  int channels, int quality = 85);

// Encode raw RGB/RGBA pixels as WebP
// Input: row-major pixel data, width * height * channels bytes
// Returns WebP-encoded bytes
std::vector<uint8_t> encode_webp(const uint8_t *data, int width, int height,
                                  int channels, int quality = 85);

// Encode row-major 8-bit pixels as PNG (zlib/deflate IDAT; no external libpng).
// channels: 1 (grayscale), 3 (RGB), or 4 (RGBA).
std::vector<uint8_t> encode_png(const uint8_t *pixels, int width, int height, int channels);

// Interleave sequential band data into BIP (Band Interleaved by Pixel) format
// Input: vector of per-band raw byte buffers (all same size = width * height * dtype_size)
// Output: interleaved bytes [R0,G0,B0,R1,G1,B1,...] where each element is dtype_size bytes
std::vector<uint8_t> interleave_bands(const std::vector<std::vector<uint8_t>> &bands,
                                       int width, int height, size_t dtype_size);

// Pack one band tile (row-major) from logical double values into raw bytes matching Raquet band layout.
// Integer types: values are rounded with std::llround then clamped to the representable range.
std::vector<uint8_t> encode_band_from_doubles(const std::vector<double> &pixels, BandDataType dtype);

} // namespace raquet
} // namespace duckdb
