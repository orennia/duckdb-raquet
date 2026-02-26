#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "band_decoder.hpp"
#include "raquet_metadata.hpp"
#include "quadbin.hpp"
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

namespace duckdb {

// ============================================================================
// WKB raster pixel type constants (PostGIS RT_PIXTYPE)
// ============================================================================

static const uint8_t WKB_PT_8BSI  = 3;   // 8-bit signed integer
static const uint8_t WKB_PT_8BUI  = 4;   // 8-bit unsigned integer
static const uint8_t WKB_PT_16BSI = 5;   // 16-bit signed integer
static const uint8_t WKB_PT_16BUI = 6;   // 16-bit unsigned integer
static const uint8_t WKB_PT_32BSI = 7;   // 32-bit signed integer
static const uint8_t WKB_PT_32BUI = 8;   // 32-bit unsigned integer
static const uint8_t WKB_PT_32BF  = 10;  // 32-bit float
static const uint8_t WKB_PT_64BF  = 11;  // 64-bit float

// Band flag: bit 5 signals that a nodata value is present
static const uint8_t WKB_BANDTYPE_FLAG_HASNODATA = (1 << 5);

// Map raquet dtype to the closest WKB pixel type
static uint8_t dtype_to_wkb_pixtype(raquet::BandDataType dtype) {
    switch (dtype) {
    case raquet::BandDataType::UINT8:   return WKB_PT_8BUI;
    case raquet::BandDataType::INT8:    return WKB_PT_8BSI;
    case raquet::BandDataType::UINT16:  return WKB_PT_16BUI;
    case raquet::BandDataType::INT16:   return WKB_PT_16BSI;
    case raquet::BandDataType::UINT32:  return WKB_PT_32BUI;
    case raquet::BandDataType::INT32:   return WKB_PT_32BSI;
    case raquet::BandDataType::UINT64:  return WKB_PT_64BF;  // no uint64 in WKB, promote to float64
    case raquet::BandDataType::INT64:   return WKB_PT_64BF;  // no int64 in WKB, promote to float64
    case raquet::BandDataType::FLOAT16: return WKB_PT_32BF;  // no float16 in WKB, promote to float32
    case raquet::BandDataType::FLOAT32: return WKB_PT_32BF;
    case raquet::BandDataType::FLOAT64: return WKB_PT_64BF;
    }
    return WKB_PT_8BUI;
}

// Returns true when the dtype needs pixel-by-pixel conversion before output
static bool dtype_needs_conversion(raquet::BandDataType dtype) {
    return dtype == raquet::BandDataType::FLOAT16 ||
           dtype == raquet::BandDataType::UINT64  ||
           dtype == raquet::BandDataType::INT64;
}

// Byte size for a WKB pixel type
static size_t wkb_pixtype_size(uint8_t pixtype) {
    switch (pixtype) {
    case WKB_PT_8BSI:
    case WKB_PT_8BUI:  return 1;
    case WKB_PT_16BSI:
    case WKB_PT_16BUI: return 2;
    case WKB_PT_32BSI:
    case WKB_PT_32BUI:
    case WKB_PT_32BF:  return 4;
    case WKB_PT_64BF:  return 8;
    default:           return 1;
    }
}

// ============================================================================
// WKB raster serialisation helpers
// ============================================================================

// Append a little-endian value to the output buffer
template<typename T>
static void wkb_write_le(std::vector<uint8_t> &buf, T val) {
    uint8_t bytes[sizeof(T)];
    memcpy(bytes, &val, sizeof(T));
    buf.insert(buf.end(), bytes, bytes + sizeof(T));
}

// Encode a nodata double value into the WKB pixel type and append it
static void wkb_write_nodata(std::vector<uint8_t> &buf, uint8_t pixtype, double nodata) {
    switch (pixtype) {
    case WKB_PT_8BUI: {
        uint8_t v = static_cast<uint8_t>(nodata);
        buf.push_back(v);
        break;
    }
    case WKB_PT_8BSI: {
        int8_t v = static_cast<int8_t>(nodata);
        buf.push_back(static_cast<uint8_t>(v));
        break;
    }
    case WKB_PT_16BUI:
        wkb_write_le<uint16_t>(buf, static_cast<uint16_t>(nodata));
        break;
    case WKB_PT_16BSI:
        wkb_write_le<int16_t>(buf, static_cast<int16_t>(nodata));
        break;
    case WKB_PT_32BUI:
        wkb_write_le<uint32_t>(buf, static_cast<uint32_t>(nodata));
        break;
    case WKB_PT_32BSI:
        wkb_write_le<int32_t>(buf, static_cast<int32_t>(nodata));
        break;
    case WKB_PT_32BF:
        wkb_write_le<float>(buf, static_cast<float>(nodata));
        break;
    case WKB_PT_64BF:
        wkb_write_le<double>(buf, nodata);
        break;
    default:
        buf.push_back(static_cast<uint8_t>(nodata));
        break;
    }
}

// ============================================================================
// Core WKB raster builder
// ============================================================================

// Build a WKB raster binary from a single decoded/raw band.
//
// WKB raster layout (little-endian):
//   Header  (61 bytes):
//     1  byte   byte order   (0x01 = LE)
//     2  bytes  version      (0x0000)
//     2  bytes  nBands       (1)
//     8  bytes  scaleX       (pixel width in degrees)
//     8  bytes  scaleY       (pixel height in degrees, negative)
//     8  bytes  ipX          (upper-left longitude)
//     8  bytes  ipY          (upper-left latitude)
//     8  bytes  skewX        (0.0)
//     8  bytes  skewY        (0.0)
//     4  bytes  srid         (int32)
//     2  bytes  width        (pixels wide)
//     2  bytes  height       (pixels tall)
//   Per band:
//     1  byte   flags        (pixtype | WKB_BANDTYPE_FLAG_HASNODATA if nodata set)
//     N  bytes  nodata       (only when has_nodata; N = wkb_pixtype_size(pixtype))
//     W*H*N bytes pixel data

static std::vector<uint8_t> build_wkb_raster(
    uint64_t block,
    int width, int height,
    raquet::BandDataType dtype,
    bool compressed,
    const uint8_t *band_raw, size_t band_size,
    bool has_nodata, double nodata,
    int32_t srid)
{
    // Geographic transform from quadbin tile
    int tile_x, tile_y, tile_z;
    quadbin::cell_to_tile(block, tile_x, tile_y, tile_z);

    double min_lon, min_lat, max_lon, max_lat;
    quadbin::tile_to_bbox_wgs84(tile_x, tile_y, tile_z,
                                 min_lon, min_lat, max_lon, max_lat);

    double scale_x =  (max_lon - min_lon) / width;
    double scale_y = -((max_lat - min_lat) / height);  // negative = north-up
    double ip_x = min_lon;
    double ip_y = max_lat;

    // WKB pixel type
    uint8_t pixtype = dtype_to_wkb_pixtype(dtype);
    size_t wkb_ps = wkb_pixtype_size(pixtype);
    size_t num_pixels = static_cast<size_t>(width) * height;

    // Decompress if needed
    std::vector<uint8_t> decompressed;
    const uint8_t *pixel_data = band_raw;
    size_t pixel_data_size = band_size;

    if (compressed) {
        decompressed = raquet::decompress_gzip(band_raw, band_size);
        pixel_data = decompressed.data();
        pixel_data_size = decompressed.size();
    }

    // Convert float16 / uint64 / int64 pixels to the promoted WKB type
    std::vector<uint8_t> converted;
    if (dtype_needs_conversion(dtype)) {
        converted.resize(num_pixels * wkb_ps);
        for (size_t i = 0; i < num_pixels; i++) {
            double val = raquet::get_pixel_value(pixel_data, pixel_data_size, i, dtype);
            if (pixtype == WKB_PT_32BF) {
                float fval = static_cast<float>(val);
                memcpy(converted.data() + i * 4, &fval, 4);
            } else {
                memcpy(converted.data() + i * 8, &val, 8);
            }
        }
        pixel_data = converted.data();
        pixel_data_size = converted.size();
    }

    size_t expected_pixel_bytes = num_pixels * wkb_ps;
    if (pixel_data_size < expected_pixel_bytes) {
        throw InvalidInputException("ST_AsWKB: band data too small for declared dimensions");
    }

    // Assemble WKB binary
    std::vector<uint8_t> wkb;
    size_t nodata_bytes = has_nodata ? wkb_ps : 0;
    wkb.reserve(61 + 1 + nodata_bytes + expected_pixel_bytes);

    // --- Header ---
    wkb.push_back(1);                                    // byte order: little endian
    wkb_write_le<uint16_t>(wkb, 0);                     // version
    wkb_write_le<uint16_t>(wkb, 1);                     // nBands = 1
    wkb_write_le<double>(wkb, scale_x);
    wkb_write_le<double>(wkb, scale_y);
    wkb_write_le<double>(wkb, ip_x);
    wkb_write_le<double>(wkb, ip_y);
    wkb_write_le<double>(wkb, 0.0);                     // skewX
    wkb_write_le<double>(wkb, 0.0);                     // skewY
    wkb_write_le<int32_t>(wkb, srid);
    wkb_write_le<uint16_t>(wkb, static_cast<uint16_t>(width));
    wkb_write_le<uint16_t>(wkb, static_cast<uint16_t>(height));

    // --- Band ---
    uint8_t flags = pixtype;
    if (has_nodata) {
        flags |= WKB_BANDTYPE_FLAG_HASNODATA;
    }
    wkb.push_back(flags);

    if (has_nodata) {
        wkb_write_nodata(wkb, pixtype, nodata);
    }

    wkb.insert(wkb.end(), pixel_data, pixel_data + expected_pixel_bytes);

    return wkb;
}

// ============================================================================
// ST_AsWKB(band, block, dtype, width, height, compression) -> BLOB
// ============================================================================

static void STAsWKBFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());
    args.data[3].Flatten(args.size());
    args.data[4].Flatten(args.size());
    args.data[5].Flatten(args.size());

    auto band_data        = FlatVector::GetData<string_t>(args.data[0]);
    auto block_data       = FlatVector::GetData<uint64_t>(args.data[1]);
    auto dtype_data       = FlatVector::GetData<string_t>(args.data[2]);
    auto width_data       = FlatVector::GetData<int32_t>(args.data[3]);
    auto height_data      = FlatVector::GetData<int32_t>(args.data[4]);
    auto compression_data = FlatVector::GetData<string_t>(args.data[5]);

    auto &band_validity = FlatVector::Validity(args.data[0]);

    for (idx_t i = 0; i < args.size(); i++) {
        if (!band_validity.RowIsValid(i)) {
            FlatVector::Validity(result).SetInvalid(i);
            continue;
        }
        try {
            auto dtype = raquet::parse_dtype(dtype_data[i].GetString());
            bool compressed = (compression_data[i].GetString() == "gzip");
            auto wkb = build_wkb_raster(
                block_data[i],
                width_data[i], height_data[i],
                dtype, compressed,
                reinterpret_cast<const uint8_t *>(band_data[i].GetData()),
                band_data[i].GetSize(),
                false, 0.0,
                4326);
            FlatVector::GetData<string_t>(result)[i] = StringVector::AddStringOrBlob(
                result, reinterpret_cast<const char *>(wkb.data()), wkb.size());
        } catch (const std::exception &) {
            FlatVector::Validity(result).SetInvalid(i);
        }
    }
}

// ============================================================================
// ST_AsWKB(band, block, dtype, width, height, compression, nodata) -> BLOB
// ============================================================================

static void STAsWKBNodataFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());
    args.data[3].Flatten(args.size());
    args.data[4].Flatten(args.size());
    args.data[5].Flatten(args.size());
    args.data[6].Flatten(args.size());

    auto band_data        = FlatVector::GetData<string_t>(args.data[0]);
    auto block_data       = FlatVector::GetData<uint64_t>(args.data[1]);
    auto dtype_data       = FlatVector::GetData<string_t>(args.data[2]);
    auto width_data       = FlatVector::GetData<int32_t>(args.data[3]);
    auto height_data      = FlatVector::GetData<int32_t>(args.data[4]);
    auto compression_data = FlatVector::GetData<string_t>(args.data[5]);
    auto nodata_data      = FlatVector::GetData<double>(args.data[6]);

    auto &band_validity = FlatVector::Validity(args.data[0]);

    for (idx_t i = 0; i < args.size(); i++) {
        if (!band_validity.RowIsValid(i)) {
            FlatVector::Validity(result).SetInvalid(i);
            continue;
        }
        try {
            auto dtype = raquet::parse_dtype(dtype_data[i].GetString());
            bool compressed = (compression_data[i].GetString() == "gzip");
            auto wkb = build_wkb_raster(
                block_data[i],
                width_data[i], height_data[i],
                dtype, compressed,
                reinterpret_cast<const uint8_t *>(band_data[i].GetData()),
                band_data[i].GetSize(),
                true, nodata_data[i],
                4326);
            FlatVector::GetData<string_t>(result)[i] = StringVector::AddStringOrBlob(
                result, reinterpret_cast<const char *>(wkb.data()), wkb.size());
        } catch (const std::exception &) {
            FlatVector::Validity(result).SetInvalid(i);
        }
    }
}

// ============================================================================
// ST_AsWKB(band, block, metadata) -> BLOB
// Extracts dtype, dimensions, compression, nodata and SRID from metadata
// ============================================================================

static void STAsWKBMetadataFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());

    auto band_data     = FlatVector::GetData<string_t>(args.data[0]);
    auto block_data    = FlatVector::GetData<uint64_t>(args.data[1]);
    auto metadata_data = FlatVector::GetData<string_t>(args.data[2]);

    auto &band_validity = FlatVector::Validity(args.data[0]);

    for (idx_t i = 0; i < args.size(); i++) {
        if (!band_validity.RowIsValid(i)) {
            FlatVector::Validity(result).SetInvalid(i);
            continue;
        }
        try {
            auto meta = raquet::parse_metadata(metadata_data[i].GetString());
            auto dtype = raquet::parse_dtype(meta.bands.empty() ? "uint8" : meta.bands[0].second);
            bool compressed = (meta.compression == "gzip");
            bool has_nodata = !meta.band_info.empty() && meta.band_info[0].has_nodata;
            double nodata   = has_nodata ? meta.band_info[0].nodata : 0.0;

            // Extract numeric SRID from "EPSG:NNNN" or default to 4326
            int32_t srid = 4326;
            if (!meta.crs.empty() && meta.crs.rfind("EPSG:", 0) == 0) {
                try { srid = std::stoi(meta.crs.substr(5)); } catch (...) {}
            }

            auto wkb = build_wkb_raster(
                block_data[i],
                meta.block_width, meta.block_height,
                dtype, compressed,
                reinterpret_cast<const uint8_t *>(band_data[i].GetData()),
                band_data[i].GetSize(),
                has_nodata, nodata,
                srid);
            FlatVector::GetData<string_t>(result)[i] = StringVector::AddStringOrBlob(
                result, reinterpret_cast<const char *>(wkb.data()), wkb.size());
        } catch (const std::exception &) {
            FlatVector::Validity(result).SetInvalid(i);
        }
    }
}

// ============================================================================
// ST_AsWKB(band, block, metadata, band_index) -> BLOB
// Multi-band variant: selects dtype and nodata by band_index
// ============================================================================

static void STAsWKBMetadataBandFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());
    args.data[3].Flatten(args.size());

    auto band_data      = FlatVector::GetData<string_t>(args.data[0]);
    auto block_data     = FlatVector::GetData<uint64_t>(args.data[1]);
    auto metadata_data  = FlatVector::GetData<string_t>(args.data[2]);
    auto band_idx_data  = FlatVector::GetData<int32_t>(args.data[3]);

    auto &band_validity = FlatVector::Validity(args.data[0]);

    for (idx_t i = 0; i < args.size(); i++) {
        if (!band_validity.RowIsValid(i)) {
            FlatVector::Validity(result).SetInvalid(i);
            continue;
        }
        try {
            auto meta = raquet::parse_metadata(metadata_data[i].GetString());
            int32_t band_idx = band_idx_data[i];
            auto dtype = raquet::parse_dtype(meta.get_band_type(band_idx));
            bool compressed = (meta.compression == "gzip");
            bool has_nodata = band_idx < static_cast<int32_t>(meta.band_info.size()) &&
                              meta.band_info[band_idx].has_nodata;
            double nodata   = has_nodata ? meta.band_info[band_idx].nodata : 0.0;

            int32_t srid = 4326;
            if (!meta.crs.empty() && meta.crs.rfind("EPSG:", 0) == 0) {
                try { srid = std::stoi(meta.crs.substr(5)); } catch (...) {}
            }

            auto wkb = build_wkb_raster(
                block_data[i],
                meta.block_width, meta.block_height,
                dtype, compressed,
                reinterpret_cast<const uint8_t *>(band_data[i].GetData()),
                band_data[i].GetSize(),
                has_nodata, nodata,
                srid);
            FlatVector::GetData<string_t>(result)[i] = StringVector::AddStringOrBlob(
                result, reinterpret_cast<const char *>(wkb.data()), wkb.size());
        } catch (const std::exception &) {
            FlatVector::Validity(result).SetInvalid(i);
        }
    }
}

// ============================================================================
// Function Registration
// ============================================================================

void RegisterAsWKBFunctions(ExtensionLoader &loader) {
    // ST_AsWKB(band, block, dtype, width, height, compression) -> BLOB
    ScalarFunction as_wkb_fn("ST_AsWKB",
        {LogicalType::BLOB, LogicalType::UBIGINT,
         LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER,
         LogicalType::VARCHAR},
        LogicalType::BLOB,
        STAsWKBFunction);
    loader.RegisterFunction(as_wkb_fn);

    // ST_AsWKB(band, block, dtype, width, height, compression, nodata) -> BLOB
    ScalarFunction as_wkb_nodata_fn("ST_AsWKB",
        {LogicalType::BLOB, LogicalType::UBIGINT,
         LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER,
         LogicalType::VARCHAR, LogicalType::DOUBLE},
        LogicalType::BLOB,
        STAsWKBNodataFunction);
    loader.RegisterFunction(as_wkb_nodata_fn);

    // ST_AsWKB(band, block, metadata) -> BLOB
    ScalarFunction as_wkb_meta_fn("ST_AsWKB",
        {LogicalType::BLOB, LogicalType::UBIGINT, LogicalType::VARCHAR},
        LogicalType::BLOB,
        STAsWKBMetadataFunction);
    loader.RegisterFunction(as_wkb_meta_fn);

    // ST_AsWKB(band, block, metadata, band_index) -> BLOB
    ScalarFunction as_wkb_meta_band_fn("ST_AsWKB",
        {LogicalType::BLOB, LogicalType::UBIGINT,
         LogicalType::VARCHAR, LogicalType::INTEGER},
        LogicalType::BLOB,
        STAsWKBMetadataBandFunction);
    loader.RegisterFunction(as_wkb_meta_band_fn);
}

} // namespace duckdb
