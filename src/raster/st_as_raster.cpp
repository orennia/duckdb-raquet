#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "band_decoder.hpp"
#include "quadbin.hpp"
#include <cmath>
#include <limits>
#include <cstring>
#include <string>
#include "yyjson.hpp"

using namespace duckdb_yyjson;

namespace duckdb {

// ============================================================================
// Geometry helper functions (for point-in-polygon tests)
// ============================================================================

// Check if a point is inside a polygon ring using ray casting
static bool AsRasterPointInRing(double px, double py, const double *ring_coords, uint32_t num_points) {
    bool inside = false;
    double x1, y1, x2, y2;

    for (uint32_t i = 0, j = num_points - 1; i < num_points; j = i++) {
        x1 = ring_coords[i * 2];
        y1 = ring_coords[i * 2 + 1];
        x2 = ring_coords[j * 2];
        y2 = ring_coords[j * 2 + 1];

        if (((y1 > py) != (y2 > py)) &&
            (px < (x2 - x1) * (py - y1) / (y2 - y1) + x1)) {
            inside = !inside;
        }
    }
    return inside;
}

// Check if a point is inside a polygon (handles holes)
static bool AsRasterPointInPolygonImpl(double px, double py, const uint8_t *data, idx_t &offset, idx_t size) {
    if (size < offset + 4) return false;
    uint32_t num_rings;
    memcpy(&num_rings, data + offset, 4);
    offset += 4;

    if (num_rings == 0) return false;

    // First ring is outer boundary
    if (size < offset + 4) return false;
    uint32_t outer_points;
    memcpy(&outer_points, data + offset, 4);
    offset += 4;

    if (size < offset + outer_points * 16) return false;

    std::vector<double> outer_coords(outer_points * 2);
    for (uint32_t i = 0; i < outer_points; i++) {
        memcpy(&outer_coords[i * 2], data + offset + i * 16, 8);
        memcpy(&outer_coords[i * 2 + 1], data + offset + i * 16 + 8, 8);
    }
    offset += outer_points * 16;

    if (!AsRasterPointInRing(px, py, outer_coords.data(), outer_points)) {
        // Skip remaining rings
        for (uint32_t r = 1; r < num_rings; r++) {
            if (size < offset + 4) return false;
            uint32_t ring_points;
            memcpy(&ring_points, data + offset, 4);
            offset += 4 + ring_points * 16;
        }
        return false;
    }

    // Check holes (if inside a hole, point is outside polygon)
    for (uint32_t r = 1; r < num_rings; r++) {
        if (size < offset + 4) return false;
        uint32_t hole_points;
        memcpy(&hole_points, data + offset, 4);
        offset += 4;

        if (size < offset + hole_points * 16) return false;

        std::vector<double> hole_coords(hole_points * 2);
        for (uint32_t i = 0; i < hole_points; i++) {
            memcpy(&hole_coords[i * 2], data + offset + i * 16, 8);
            memcpy(&hole_coords[i * 2 + 1], data + offset + i * 16 + 8, 8);
        }
        offset += hole_points * 16;

        if (AsRasterPointInRing(px, py, hole_coords.data(), hole_points)) {
            return false; // Inside a hole
        }
    }

    return true;
}

// Check if a point is inside a geometry (polygon or multipolygon)
static bool AsRasterPointInGeometry(double px, double py, const string_t &geom) {
    const uint8_t *data = reinterpret_cast<const uint8_t *>(geom.GetData());
    idx_t size = geom.GetSize();

    if (size < 5) return false;

    uint32_t geom_type;
    memcpy(&geom_type, data + 1, 4);

    uint32_t base_type = geom_type & 0xFF;
    idx_t offset = 5;

    // Handle SRID prefix
    if (geom_type & 0x20000000) {
        offset += 4;
    }

    if (base_type == 3) { // POLYGON
        return AsRasterPointInPolygonImpl(px, py, data, offset, size);
    } else if (base_type == 6) { // MULTIPOLYGON
        if (size < offset + 4) return false;
        uint32_t num_polygons;
        memcpy(&num_polygons, data + offset, 4);
        offset += 4;

        for (uint32_t poly = 0; poly < num_polygons; poly++) {
            if (size < offset + 5) return false;
            offset += 5;

            uint32_t poly_type;
            memcpy(&poly_type, data + offset - 4, 4);
            if (poly_type & 0x20000000) {
                offset += 4;
            }

            idx_t poly_offset = offset;
            if (AsRasterPointInPolygonImpl(px, py, data, poly_offset, size)) {
                return true;
            }
            offset = poly_offset;
        }
        return false;
    }

    return false;
}

// Extract bounding box from geometry
static bool AsRasterExtractGeometryBBox(const string_t &geom, double &min_x, double &min_y,
                                         double &max_x, double &max_y) {
    const uint8_t *data = reinterpret_cast<const uint8_t *>(geom.GetData());
    idx_t size = geom.GetSize();

    if (size < 5) return false;

    uint32_t geom_type;
    memcpy(&geom_type, data + 1, 4);

    uint32_t base_type = geom_type & 0xFF;
    idx_t offset = 5;

    if (geom_type & 0x20000000) {
        offset += 4;
    }

    min_x = std::numeric_limits<double>::max();
    min_y = std::numeric_limits<double>::max();
    max_x = std::numeric_limits<double>::lowest();
    max_y = std::numeric_limits<double>::lowest();

    auto update_bounds = [&](double x, double y) {
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    };

    if (base_type == 3) { // POLYGON
        if (size < offset + 4) return false;
        uint32_t num_rings;
        memcpy(&num_rings, data + offset, 4);
        offset += 4;

        for (uint32_t r = 0; r < num_rings; r++) {
            if (size < offset + 4) return false;
            uint32_t num_points;
            memcpy(&num_points, data + offset, 4);
            offset += 4;

            if (size < offset + num_points * 16) return false;
            for (uint32_t p = 0; p < num_points; p++) {
                double x, y;
                memcpy(&x, data + offset + p * 16, 8);
                memcpy(&y, data + offset + p * 16 + 8, 8);
                update_bounds(x, y);
            }
            offset += num_points * 16;
        }
        return true;
    } else if (base_type == 6) { // MULTIPOLYGON
        if (size < offset + 4) return false;
        uint32_t num_polygons;
        memcpy(&num_polygons, data + offset, 4);
        offset += 4;

        for (uint32_t poly = 0; poly < num_polygons; poly++) {
            if (size < offset + 5) return false;
            offset += 5;

            uint32_t poly_type;
            memcpy(&poly_type, data + offset - 4, 4);
            if (poly_type & 0x20000000) {
                offset += 4;
            }

            if (size < offset + 4) return false;
            uint32_t num_rings;
            memcpy(&num_rings, data + offset, 4);
            offset += 4;

            for (uint32_t r = 0; r < num_rings; r++) {
                if (size < offset + 4) return false;
                uint32_t num_points;
                memcpy(&num_points, data + offset, 4);
                offset += 4;

                if (size < offset + num_points * 16) return false;
                for (uint32_t p = 0; p < num_points; p++) {
                    double x, y;
                    memcpy(&x, data + offset + p * 16, 8);
                    memcpy(&y, data + offset + p * 16 + 8, 8);
                    update_bounds(x, y);
                }
                offset += num_points * 16;
            }
        }
        return true;
    }

    return false;
}

// ============================================================================
// Pixel value encoding
// ============================================================================

// Convert double to IEEE 754 half-precision float16
static uint16_t double_to_float16(double val) {
    if (std::isnan(val)) return 0x7E00;
    if (std::isinf(val)) return (val > 0) ? 0x7C00 : 0xFC00;

    float f = static_cast<float>(val);
    uint32_t bits;
    memcpy(&bits, &f, 4);

    uint16_t sign = static_cast<uint16_t>((bits >> 16) & 0x8000);
    int32_t exp = static_cast<int32_t>((bits >> 23) & 0xFF) - 127 + 15;
    uint32_t mant = bits & 0x7FFFFF;

    if (exp <= 0) {
        if (exp < -10) return sign; // Too small, flush to zero
        mant = (mant | 0x800000) >> (1 - exp);
        return sign | static_cast<uint16_t>(mant >> 13);
    }
    if (exp >= 31) return sign | 0x7C00; // Overflow to infinity

    return sign | static_cast<uint16_t>(exp << 10) | static_cast<uint16_t>(mant >> 13);
}

// Write a typed pixel value into the output buffer at pixel_index
static void set_pixel_value(uint8_t *buffer, size_t data_size, size_t pixel_index,
                             raquet::BandDataType dtype, double value) {
    size_t elem_size = raquet::dtype_size(dtype);
    size_t byte_offset = pixel_index * elem_size;

    if (byte_offset + elem_size > data_size) return;

    switch (dtype) {
    case raquet::BandDataType::UINT8: {
        uint8_t v = (value < 0.0) ? 0 : (value > 255.0) ? 255 : static_cast<uint8_t>(value);
        buffer[byte_offset] = v;
        break;
    }
    case raquet::BandDataType::INT8: {
        int8_t v = (value < -128.0) ? -128 : (value > 127.0) ? 127 : static_cast<int8_t>(value);
        memcpy(buffer + byte_offset, &v, 1);
        break;
    }
    case raquet::BandDataType::UINT16: {
        uint16_t v = (value < 0.0) ? 0 : (value > 65535.0) ? 65535 : static_cast<uint16_t>(value);
        memcpy(buffer + byte_offset, &v, 2);
        break;
    }
    case raquet::BandDataType::INT16: {
        int16_t v = (value < -32768.0) ? -32768 : (value > 32767.0) ? 32767 : static_cast<int16_t>(value);
        memcpy(buffer + byte_offset, &v, 2);
        break;
    }
    case raquet::BandDataType::UINT32: {
        uint32_t v = (value < 0.0) ? 0 : static_cast<uint32_t>(value);
        memcpy(buffer + byte_offset, &v, 4);
        break;
    }
    case raquet::BandDataType::INT32: {
        int32_t v = static_cast<int32_t>(value);
        memcpy(buffer + byte_offset, &v, 4);
        break;
    }
    case raquet::BandDataType::UINT64: {
        uint64_t v = (value < 0.0) ? 0 : static_cast<uint64_t>(value);
        memcpy(buffer + byte_offset, &v, 8);
        break;
    }
    case raquet::BandDataType::INT64: {
        int64_t v = static_cast<int64_t>(value);
        memcpy(buffer + byte_offset, &v, 8);
        break;
    }
    case raquet::BandDataType::FLOAT16: {
        uint16_t v = double_to_float16(value);
        memcpy(buffer + byte_offset, &v, 2);
        break;
    }
    case raquet::BandDataType::FLOAT32: {
        float v = static_cast<float>(value);
        memcpy(buffer + byte_offset, &v, 4);
        break;
    }
    case raquet::BandDataType::FLOAT64: {
        memcpy(buffer + byte_offset, &value, 8);
        break;
    }
    }
}

// ============================================================================
// Metadata parsing
// ============================================================================

struct AsRasterMetadata {
    std::string dtype;
    int block_width;
    int block_height;

    static AsRasterMetadata Parse(const std::string &json_str) {
        AsRasterMetadata meta;
        meta.dtype = "uint8";
        meta.block_width = 256;
        meta.block_height = 256;

        yyjson_doc *doc = yyjson_read(json_str.c_str(), json_str.length(), 0);
        if (!doc) {
            throw InvalidInputException("ST_AsRaster: Invalid metadata JSON");
        }

        yyjson_val *root = yyjson_doc_get_root(doc);
        if (!yyjson_is_obj(root)) {
            yyjson_doc_free(doc);
            throw InvalidInputException("ST_AsRaster: Metadata must be a JSON object");
        }

        // Parse tiling dimensions
        yyjson_val *tiling_val = yyjson_obj_get(root, "tiling");
        if (tiling_val && yyjson_is_obj(tiling_val)) {
            yyjson_val *width_val = yyjson_obj_get(tiling_val, "block_width");
            if (width_val && yyjson_is_int(width_val)) {
                meta.block_width = yyjson_get_int(width_val);
            }
            yyjson_val *height_val = yyjson_obj_get(tiling_val, "block_height");
            if (height_val && yyjson_is_int(height_val)) {
                meta.block_height = yyjson_get_int(height_val);
            }
        }

        // Parse band data type from first band
        yyjson_val *bands_val = yyjson_obj_get(root, "bands");
        if (bands_val && yyjson_is_arr(bands_val)) {
            size_t idx, max;
            yyjson_val *band;
            yyjson_arr_foreach(bands_val, idx, max, band) {
                if (yyjson_is_obj(band)) {
                    yyjson_val *type_val = yyjson_obj_get(band, "type");
                    if (type_val && yyjson_is_str(type_val)) {
                        meta.dtype = yyjson_get_str(type_val);
                    }
                }
                break; // Only need the first band's type
            }
        }

        yyjson_doc_free(doc);
        return meta;
    }
};

// ============================================================================
// Core rasterization logic
// ============================================================================

// Rasterize a geometry within the bounds of a quadbin tile
// Returns an uncompressed BLOB of pixel values (row-major order)
static std::vector<uint8_t> rasterize_geometry(const string_t &geom, uint64_t block,
                                                 int width, int height,
                                                 raquet::BandDataType dtype,
                                                 double value, double nodata) {
    // Get tile geographic bounds
    int tile_x, tile_y, tile_z;
    quadbin::cell_to_tile(block, tile_x, tile_y, tile_z);

    double tile_min_lon, tile_min_lat, tile_max_lon, tile_max_lat;
    quadbin::tile_to_bbox_wgs84(tile_x, tile_y, tile_z,
                                 tile_min_lon, tile_min_lat, tile_max_lon, tile_max_lat);

    size_t pixel_count = static_cast<size_t>(width) * height;
    size_t elem_size = raquet::dtype_size(dtype);
    size_t data_size = pixel_count * elem_size;

    std::vector<uint8_t> result(data_size);

    // Initialize all pixels to nodata
    for (size_t i = 0; i < pixel_count; i++) {
        set_pixel_value(result.data(), data_size, i, dtype, nodata);
    }

    // Get geometry bounding box for fast rejection
    double geom_min_x, geom_min_y, geom_max_x, geom_max_y;
    bool has_bbox = AsRasterExtractGeometryBBox(geom, geom_min_x, geom_min_y, geom_max_x, geom_max_y);

    // Fast reject: geometry bbox doesn't overlap tile
    if (has_bbox && (tile_max_lon < geom_min_x || tile_min_lon > geom_max_x ||
                     tile_max_lat < geom_min_y || tile_min_lat > geom_max_y)) {
        return result; // All nodata
    }

    double pixel_width = (tile_max_lon - tile_min_lon) / width;
    double pixel_height = (tile_max_lat - tile_min_lat) / height;

    for (int py = 0; py < height; py++) {
        for (int px = 0; px < width; px++) {
            // Pixel center coordinates (lon/lat)
            double pixel_lon = tile_min_lon + (px + 0.5) * pixel_width;
            double pixel_lat = tile_max_lat - (py + 0.5) * pixel_height;

            // Quick bbox check
            if (has_bbox && (pixel_lon < geom_min_x || pixel_lon > geom_max_x ||
                             pixel_lat < geom_min_y || pixel_lat > geom_max_y)) {
                continue;
            }

            if (AsRasterPointInGeometry(pixel_lon, pixel_lat, geom)) {
                size_t pixel_index = static_cast<size_t>(py) * width + px;
                set_pixel_value(result.data(), data_size, pixel_index, dtype, value);
            }
        }
    }

    return result;
}

// ============================================================================
// ST_AsRaster(geometry, block, width, height, dtype) -> BLOB
// Rasterizes geometry within tile; inside pixels = 1, outside = 0
// ============================================================================

static void STAsRasterFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());
    args.data[3].Flatten(args.size());
    args.data[4].Flatten(args.size());

    auto geom_data = FlatVector::GetData<string_t>(args.data[0]);
    auto block_data = FlatVector::GetData<uint64_t>(args.data[1]);
    auto width_data = FlatVector::GetData<int32_t>(args.data[2]);
    auto height_data = FlatVector::GetData<int32_t>(args.data[3]);
    auto dtype_data = FlatVector::GetData<string_t>(args.data[4]);

    auto &geom_validity = FlatVector::Validity(args.data[0]);

    for (idx_t i = 0; i < args.size(); i++) {
        if (!geom_validity.RowIsValid(i)) {
            FlatVector::Validity(result).SetInvalid(i);
            continue;
        }

        try {
            auto dtype = raquet::parse_dtype(dtype_data[i].GetString());
            auto pixels = rasterize_geometry(geom_data[i], block_data[i],
                                              width_data[i], height_data[i],
                                              dtype, 1.0, 0.0);
            FlatVector::GetData<string_t>(result)[i] = StringVector::AddStringOrBlob(result,
                reinterpret_cast<const char *>(pixels.data()), pixels.size());
        } catch (const std::exception &) {
            FlatVector::Validity(result).SetInvalid(i);
        }
    }
}

// ============================================================================
// ST_AsRaster(geometry, block, width, height, dtype, value, nodata) -> BLOB
// Rasterizes geometry with explicit pixel value and nodata
// ============================================================================

static void STAsRasterValueFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());
    args.data[3].Flatten(args.size());
    args.data[4].Flatten(args.size());
    args.data[5].Flatten(args.size());
    args.data[6].Flatten(args.size());

    auto geom_data = FlatVector::GetData<string_t>(args.data[0]);
    auto block_data = FlatVector::GetData<uint64_t>(args.data[1]);
    auto width_data = FlatVector::GetData<int32_t>(args.data[2]);
    auto height_data = FlatVector::GetData<int32_t>(args.data[3]);
    auto dtype_data = FlatVector::GetData<string_t>(args.data[4]);
    auto value_data = FlatVector::GetData<double>(args.data[5]);
    auto nodata_data = FlatVector::GetData<double>(args.data[6]);

    auto &geom_validity = FlatVector::Validity(args.data[0]);

    for (idx_t i = 0; i < args.size(); i++) {
        if (!geom_validity.RowIsValid(i)) {
            FlatVector::Validity(result).SetInvalid(i);
            continue;
        }

        try {
            auto dtype = raquet::parse_dtype(dtype_data[i].GetString());
            auto pixels = rasterize_geometry(geom_data[i], block_data[i],
                                              width_data[i], height_data[i],
                                              dtype, value_data[i], nodata_data[i]);
            FlatVector::GetData<string_t>(result)[i] = StringVector::AddStringOrBlob(result,
                reinterpret_cast<const char *>(pixels.data()), pixels.size());
        } catch (const std::exception &) {
            FlatVector::Validity(result).SetInvalid(i);
        }
    }
}

// ============================================================================
// ST_AsRaster(geometry, block, metadata) -> BLOB
// Rasterizes geometry using tile dimensions and dtype from Raquet metadata
// Inside pixels = 1, outside = 0
// ============================================================================

static void STAsRasterMetadataFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());

    auto geom_data = FlatVector::GetData<string_t>(args.data[0]);
    auto block_data = FlatVector::GetData<uint64_t>(args.data[1]);
    auto metadata_data = FlatVector::GetData<string_t>(args.data[2]);

    auto &geom_validity = FlatVector::Validity(args.data[0]);

    for (idx_t i = 0; i < args.size(); i++) {
        if (!geom_validity.RowIsValid(i)) {
            FlatVector::Validity(result).SetInvalid(i);
            continue;
        }

        try {
            auto meta = AsRasterMetadata::Parse(metadata_data[i].GetString());
            auto dtype = raquet::parse_dtype(meta.dtype);
            auto pixels = rasterize_geometry(geom_data[i], block_data[i],
                                              meta.block_width, meta.block_height,
                                              dtype, 1.0, 0.0);
            FlatVector::GetData<string_t>(result)[i] = StringVector::AddStringOrBlob(result,
                reinterpret_cast<const char *>(pixels.data()), pixels.size());
        } catch (const std::exception &) {
            FlatVector::Validity(result).SetInvalid(i);
        }
    }
}

// ============================================================================
// ST_AsRaster(geometry, block, metadata, value, nodata) -> BLOB
// Rasterizes geometry using metadata dimensions/dtype with explicit pixel values
// ============================================================================

static void STAsRasterMetadataValueFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());
    args.data[3].Flatten(args.size());
    args.data[4].Flatten(args.size());

    auto geom_data = FlatVector::GetData<string_t>(args.data[0]);
    auto block_data = FlatVector::GetData<uint64_t>(args.data[1]);
    auto metadata_data = FlatVector::GetData<string_t>(args.data[2]);
    auto value_data = FlatVector::GetData<double>(args.data[3]);
    auto nodata_data = FlatVector::GetData<double>(args.data[4]);

    auto &geom_validity = FlatVector::Validity(args.data[0]);

    for (idx_t i = 0; i < args.size(); i++) {
        if (!geom_validity.RowIsValid(i)) {
            FlatVector::Validity(result).SetInvalid(i);
            continue;
        }

        try {
            auto meta = AsRasterMetadata::Parse(metadata_data[i].GetString());
            auto dtype = raquet::parse_dtype(meta.dtype);
            auto pixels = rasterize_geometry(geom_data[i], block_data[i],
                                              meta.block_width, meta.block_height,
                                              dtype, value_data[i], nodata_data[i]);
            FlatVector::GetData<string_t>(result)[i] = StringVector::AddStringOrBlob(result,
                reinterpret_cast<const char *>(pixels.data()), pixels.size());
        } catch (const std::exception &) {
            FlatVector::Validity(result).SetInvalid(i);
        }
    }
}

// ============================================================================
// Function Registration
// ============================================================================

void RegisterAsRasterFunctions(ExtensionLoader &loader) {
    // ST_AsRaster(geometry, block, width, height, dtype) -> BLOB
    // Rasterizes geometry within tile; pixels inside = 1, outside = 0
    ScalarFunction as_raster_fn("ST_AsRaster",
        {LogicalType::GEOMETRY(), LogicalType::UBIGINT,
         LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::VARCHAR},
        LogicalType::BLOB,
        STAsRasterFunction);
    loader.RegisterFunction(as_raster_fn);

    // ST_AsRaster(geometry, block, width, height, dtype, value, nodata) -> BLOB
    // Rasterizes geometry with explicit pixel value and nodata value
    ScalarFunction as_raster_value_fn("ST_AsRaster",
        {LogicalType::GEOMETRY(), LogicalType::UBIGINT,
         LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::VARCHAR,
         LogicalType::DOUBLE, LogicalType::DOUBLE},
        LogicalType::BLOB,
        STAsRasterValueFunction);
    loader.RegisterFunction(as_raster_value_fn);

    // ST_AsRaster(geometry, block, metadata) -> BLOB
    // Rasterizes geometry using Raquet metadata for tile dimensions and dtype
    ScalarFunction as_raster_meta_fn("ST_AsRaster",
        {LogicalType::GEOMETRY(), LogicalType::UBIGINT, LogicalType::VARCHAR},
        LogicalType::BLOB,
        STAsRasterMetadataFunction);
    loader.RegisterFunction(as_raster_meta_fn);

    // ST_AsRaster(geometry, block, metadata, value, nodata) -> BLOB
    // Rasterizes geometry using metadata dimensions/dtype with explicit pixel values
    ScalarFunction as_raster_meta_value_fn("ST_AsRaster",
        {LogicalType::GEOMETRY(), LogicalType::UBIGINT, LogicalType::VARCHAR,
         LogicalType::DOUBLE, LogicalType::DOUBLE},
        LogicalType::BLOB,
        STAsRasterMetadataValueFunction);
    loader.RegisterFunction(as_raster_meta_value_fn);
}

} // namespace duckdb
