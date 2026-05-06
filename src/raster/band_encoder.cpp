#include "band_encoder.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <zlib.h>

#ifdef RAQUET_HAS_JPEG
#include <jpeglib.h>
#endif

#ifdef RAQUET_HAS_WEBP
#include <webp/encode.h>
#endif

namespace duckdb {
namespace raquet {

// IEEE float32 → float16 (half); sufficient for round-tripping masked tiles.
static uint16_t float_to_half_bits(float f) {
	uint32_t x;
	std::memcpy(&x, &f, sizeof(x));
	uint32_t sign = (x >> 16) & 0x8000u;
	uint32_t mant = x & 0x007fffffu;
	int32_t exp = static_cast<int32_t>((x >> 23) & 0xffu);

	if (exp == 255) {
		return static_cast<uint16_t>(sign | (mant != 0 ? 0x7e00u : 0x7c00u));
	}
	if (exp == 0 && mant == 0) {
		return static_cast<uint16_t>(sign);
	}
	exp = exp - 127 + 15;
	if (exp >= 31) {
		return static_cast<uint16_t>(sign | 0x7c00u);
	}
	if (exp <= 0) {
		if (exp < -10) {
			return static_cast<uint16_t>(sign);
		}
		mant = (mant | 0x00800000u) >> static_cast<uint32_t>(1 - exp);
		return static_cast<uint16_t>(sign | (mant >> 13));
	}
	return static_cast<uint16_t>(sign | (static_cast<uint32_t>(exp) << 10) | (static_cast<uint32_t>(mant) >> 13));
}

template <typename T>
static T clamp_round_integer(double v) {
	if (std::isnan(v) || std::isinf(v)) {
		return T {0};
	}
	long double x = std::llround(static_cast<long double>(v));
	const long double lo = static_cast<long double>(std::numeric_limits<T>::min());
	const long double hi = static_cast<long double>(std::numeric_limits<T>::max());
	if (x < lo) {
		x = lo;
	}
	if (x > hi) {
		x = hi;
	}
	return static_cast<T>(x);
}

std::vector<uint8_t> encode_band_from_doubles(const std::vector<double> &pixels, BandDataType dtype) {
	const size_t elem_size = dtype_size(dtype);
	std::vector<uint8_t> out(pixels.size() * elem_size);
	for (size_t i = 0; i < pixels.size(); i++) {
		const double v = pixels[i];
		uint8_t *ptr = out.data() + i * elem_size;
		switch (dtype) {
		case BandDataType::UINT8: {
			uint8_t u = clamp_round_integer<uint8_t>(v);
			std::memcpy(ptr, &u, 1);
			break;
		}
		case BandDataType::INT8: {
			int8_t s = clamp_round_integer<int8_t>(v);
			std::memcpy(ptr, &s, 1);
			break;
		}
		case BandDataType::UINT16: {
			uint16_t u = clamp_round_integer<uint16_t>(v);
			std::memcpy(ptr, &u, 2);
			break;
		}
		case BandDataType::INT16: {
			int16_t s = clamp_round_integer<int16_t>(v);
			std::memcpy(ptr, &s, 2);
			break;
		}
		case BandDataType::UINT32: {
			uint32_t u = clamp_round_integer<uint32_t>(v);
			std::memcpy(ptr, &u, 4);
			break;
		}
		case BandDataType::INT32: {
			int32_t s = clamp_round_integer<int32_t>(v);
			std::memcpy(ptr, &s, 4);
			break;
		}
		case BandDataType::UINT64: {
			uint64_t u = 0;
			if (std::isfinite(v) && v > 0) {
				if (v >= static_cast<double>(std::numeric_limits<uint64_t>::max())) {
					u = std::numeric_limits<uint64_t>::max();
				} else {
					u = static_cast<uint64_t>(v);
				}
			}
			std::memcpy(ptr, &u, 8);
			break;
		}
		case BandDataType::INT64: {
			int64_t s = 0;
			if (std::isfinite(v)) {
				if (v >= static_cast<double>(std::numeric_limits<int64_t>::max())) {
					s = std::numeric_limits<int64_t>::max();
				} else if (v <= static_cast<double>(std::numeric_limits<int64_t>::min())) {
					s = std::numeric_limits<int64_t>::min();
				} else {
					s = static_cast<int64_t>(std::llround(v));
				}
			}
			std::memcpy(ptr, &s, 8);
			break;
		}
		case BandDataType::FLOAT16: {
			const float fv = static_cast<float>(v);
			uint16_t h = float_to_half_bits(fv);
			std::memcpy(ptr, &h, 2);
			break;
		}
		case BandDataType::FLOAT32: {
			const float fv = static_cast<float>(v);
			std::memcpy(ptr, &fv, 4);
			break;
		}
		case BandDataType::FLOAT64:
			std::memcpy(ptr, &v, 8);
			break;
		}
	}
	return out;
}

static uint32_t crc32_png(const uint8_t *buf, size_t len) {
	static uint32_t table[256];
	static bool table_init = false;
	if (!table_init) {
		for (uint32_t n = 0; n < 256; n++) {
			uint32_t c = n;
			for (int k = 0; k < 8; k++) {
				c = (c & 1U) ? (0xEDB88320U ^ (c >> 1)) : (c >> 1);
			}
			table[n] = c;
		}
		table_init = true;
	}
	uint32_t c = 0xFFFFFFFFU;
	for (size_t i = 0; i < len; i++) {
		c = table[(c ^ buf[i]) & 0xFFU] ^ (c >> 8);
	}
	return c ^ 0xFFFFFFFFU;
}

static void push_be32(std::vector<uint8_t> &v, uint32_t x) {
	v.push_back(static_cast<uint8_t>((x >> 24) & 0xFFU));
	v.push_back(static_cast<uint8_t>((x >> 16) & 0xFFU));
	v.push_back(static_cast<uint8_t>((x >> 8) & 0xFFU));
	v.push_back(static_cast<uint8_t>(x & 0xFFU));
}

static void append_png_chunk(std::vector<uint8_t> &out, const char type[4], const uint8_t *data, size_t len) {
	push_be32(out, static_cast<uint32_t>(len));
	for (int i = 0; i < 4; i++) {
		out.push_back(static_cast<uint8_t>(type[i]));
	}
	if (len > 0 && data != nullptr) {
		out.insert(out.end(), data, data + len);
	}
	std::vector<uint8_t> crc_in(4 + len);
	std::memcpy(crc_in.data(), type, 4);
	if (len > 0 && data != nullptr) {
		std::memcpy(crc_in.data() + 4, data, len);
	}
	push_be32(out, crc32_png(crc_in.data(), crc_in.size()));
}

std::vector<uint8_t> compress_gzip(const uint8_t *data, size_t size) {
    // Estimate compressed size (zlib recommends compressBound)
    uLongf compressed_size = compressBound(static_cast<uLong>(size));
    std::vector<uint8_t> compressed(compressed_size);

    int ret = compress2(compressed.data(), &compressed_size,
                        data, static_cast<uLong>(size), Z_DEFAULT_COMPRESSION);
    if (ret != Z_OK) {
        throw std::runtime_error("gzip compression failed with error code " + std::to_string(ret));
    }

    compressed.resize(compressed_size);
    return compressed;
}

std::vector<uint8_t> encode_png(const uint8_t *pixels, int width, int height, int channels) {
	if (channels != 1 && channels != 3 && channels != 4) {
		throw std::invalid_argument("encode_png: channels must be 1, 3, or 4 (got " + std::to_string(channels) + ")");
	}
	if (width <= 0 || height <= 0) {
		throw std::invalid_argument("encode_png: width and height must be positive");
	}
	const uint8_t color_type =
	    channels == 1 ? static_cast<uint8_t>(0) : (channels == 3 ? static_cast<uint8_t>(2) : static_cast<uint8_t>(6));
	const size_t bpp = static_cast<size_t>(channels);
	const size_t row_bytes = static_cast<size_t>(width) * bpp;
	std::vector<uint8_t> scanlines;
	scanlines.reserve(static_cast<size_t>(height) * (1 + row_bytes));
	for (int y = 0; y < height; y++) {
		scanlines.push_back(0); // filter type None
		const size_t row_off = static_cast<size_t>(y) * row_bytes;
		scanlines.insert(scanlines.end(), pixels + row_off, pixels + row_off + row_bytes);
	}
	// zlib (RFC 1950) stream — same compressor as Raquet tile blobs; PNG IDAT expects zlib.
	std::vector<uint8_t> idat_payload = compress_gzip(scanlines.data(), scanlines.size());

	static const uint8_t kPngSig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
	std::vector<uint8_t> png;
	png.reserve(8 + 12 + 13 + 12 + idat_payload.size() + 12);
	png.insert(png.end(), kPngSig, kPngSig + 8);

	uint8_t ihdr[13];
	auto wr32 = [](uint8_t *p, uint32_t v) {
		p[0] = static_cast<uint8_t>((v >> 24) & 0xFFU);
		p[1] = static_cast<uint8_t>((v >> 16) & 0xFFU);
		p[2] = static_cast<uint8_t>((v >> 8) & 0xFFU);
		p[3] = static_cast<uint8_t>(v & 0xFFU);
	};
	wr32(ihdr + 0, static_cast<uint32_t>(width));
	wr32(ihdr + 4, static_cast<uint32_t>(height));
	ihdr[8] = 8; // bit depth
	ihdr[9] = color_type;
	ihdr[10] = 0; // compression
	ihdr[11] = 0; // filter
	ihdr[12] = 0; // interlace
	append_png_chunk(png, "IHDR", ihdr, 13);
	append_png_chunk(png, "IDAT", idat_payload.data(), idat_payload.size());
	append_png_chunk(png, "IEND", nullptr, 0);
	return png;
}

std::vector<uint8_t> encode_jpeg(const uint8_t *data, int width, int height,
                                  int channels, int quality) {
#ifdef RAQUET_HAS_JPEG
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    // Write to memory buffer
    unsigned char *outbuffer = nullptr;
    unsigned long outsize = 0;
    jpeg_mem_dest(&cinfo, &outbuffer, &outsize);

    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = channels;

    if (channels == 1) {
        cinfo.in_color_space = JCS_GRAYSCALE;
    } else if (channels == 3) {
        cinfo.in_color_space = JCS_RGB;
    } else {
        jpeg_destroy_compress(&cinfo);
        throw std::invalid_argument("JPEG encoding supports 1 or 3 channels, got " +
                                     std::to_string(channels));
    }

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    int row_stride = width * channels;
    while (cinfo.next_scanline < cinfo.image_height) {
        const uint8_t *row_ptr = data + cinfo.next_scanline * row_stride;
        JSAMPROW row = const_cast<JSAMPROW>(row_ptr);
        jpeg_write_scanlines(&cinfo, &row, 1);
    }

    jpeg_finish_compress(&cinfo);

    std::vector<uint8_t> result(outbuffer, outbuffer + outsize);

    jpeg_destroy_compress(&cinfo);
    free(outbuffer);

    return result;
#else
    throw std::runtime_error("JPEG encoding not available (libjpeg not linked)");
#endif
}

std::vector<uint8_t> encode_webp(const uint8_t *data, int width, int height,
                                  int channels, int quality) {
#ifdef RAQUET_HAS_WEBP
    uint8_t *output = nullptr;
    size_t output_size = 0;

    if (channels == 3) {
        output_size = WebPEncodeRGB(data, width, height, width * 3,
                                     static_cast<float>(quality), &output);
    } else if (channels == 4) {
        output_size = WebPEncodeRGBA(data, width, height, width * 4,
                                      static_cast<float>(quality), &output);
    } else {
        throw std::invalid_argument("WebP encoding supports 3 or 4 channels, got " +
                                     std::to_string(channels));
    }

    if (output_size == 0 || output == nullptr) {
        if (output) WebPFree(output);
        throw std::runtime_error("WebP encoding failed");
    }

    std::vector<uint8_t> result(output, output + output_size);
    WebPFree(output);
    return result;
#else
    throw std::runtime_error("WebP encoding not available (libwebp not linked)");
#endif
}

std::vector<uint8_t> interleave_bands(const std::vector<std::vector<uint8_t>> &bands,
                                       int width, int height, size_t dtype_size) {
    size_t num_bands = bands.size();
    size_t num_pixels = static_cast<size_t>(width) * height;
    size_t total_size = num_pixels * num_bands * dtype_size;

    std::vector<uint8_t> interleaved(total_size);

    for (size_t pixel = 0; pixel < num_pixels; pixel++) {
        for (size_t band = 0; band < num_bands; band++) {
            size_t src_offset = pixel * dtype_size;
            size_t dst_offset = (pixel * num_bands + band) * dtype_size;
            std::memcpy(interleaved.data() + dst_offset,
                        bands[band].data() + src_offset,
                        dtype_size);
        }
    }

    return interleaved;
}

} // namespace raquet
} // namespace duckdb
