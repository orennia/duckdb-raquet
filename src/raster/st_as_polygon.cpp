#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "band_decoder.hpp"
#include "quadbin.hpp"
#include "raquet_metadata.hpp"
#include <cstring>
#include <vector>

namespace duckdb {

// Build a WKB polygon from bounding box coordinates (little-endian)
static std::vector<uint8_t> MakeWKBPolygon(double min_lon, double min_lat,
                                            double max_lon, double max_lat) {
    std::vector<uint8_t> wkb;
    wkb.reserve(93); // 1 + 4 + 4 + 4 + 5*16 = 93 bytes

    // Byte order: little-endian
    wkb.push_back(0x01);

    // Geometry type: Polygon (3)
    uint32_t geom_type = 3;
    wkb.insert(wkb.end(), reinterpret_cast<uint8_t*>(&geom_type),
               reinterpret_cast<uint8_t*>(&geom_type) + 4);

    // Number of rings: 1
    uint32_t num_rings = 1;
    wkb.insert(wkb.end(), reinterpret_cast<uint8_t*>(&num_rings),
               reinterpret_cast<uint8_t*>(&num_rings) + 4);

    // Number of points in ring: 5 (closed polygon)
    uint32_t num_points = 5;
    wkb.insert(wkb.end(), reinterpret_cast<uint8_t*>(&num_points),
               reinterpret_cast<uint8_t*>(&num_points) + 4);

    // Points: SW, SE, NE, NW, SW (closed ring)
    auto add_point = [&wkb](double x, double y) {
        wkb.insert(wkb.end(), reinterpret_cast<uint8_t*>(&x),
                   reinterpret_cast<uint8_t*>(&x) + 8);
        wkb.insert(wkb.end(), reinterpret_cast<uint8_t*>(&y),
                   reinterpret_cast<uint8_t*>(&y) + 8);
    };

    add_point(min_lon, min_lat); // SW
    add_point(max_lon, min_lat); // SE
    add_point(max_lon, max_lat); // NE
    add_point(min_lon, max_lat); // NW
    add_point(min_lon, min_lat); // SW (close ring)

    return wkb;
}

// ============================================================================
// ST_AsPolygon(block, metadata) -> GEOMETRY
// Returns the WGS84 bounding polygon of the raster tile (0-indexed bands).
// Similar to PostGIS: geometry ST_Polygon(raster rast, integer band_num=1)
// ============================================================================

// ST_AsPolygon(block UBIGINT, metadata VARCHAR) -> GEOMETRY
static void STAsPolygonFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());

    auto block_data = FlatVector::GetData<uint64_t>(args.data[0]);
    auto &result_mask = FlatVector::Validity(result);

    for (idx_t i = 0; i < args.size(); i++) {
        auto block = block_data[i];

        try {
            int x, y, z;
            quadbin::cell_to_tile(block, x, y, z);

            double min_lon, min_lat, max_lon, max_lat;
            quadbin::tile_to_bbox_wgs84(x, y, z, min_lon, min_lat, max_lon, max_lat);

            auto wkb = MakeWKBPolygon(min_lon, min_lat, max_lon, max_lat);
            FlatVector::GetData<string_t>(result)[i] =
                StringVector::AddStringOrBlob(result,
                    reinterpret_cast<const char*>(wkb.data()), wkb.size());
        } catch (const std::exception &e) {
            result_mask.SetInvalid(i);
        }
    }
}

// ST_AsPolygon(block UBIGINT, band BLOB, metadata VARCHAR) -> GEOMETRY
// Band data is accepted for API compatibility; band_num defaults to the first band (index 0).
static void STAsPolygonWithBandFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());

    auto block_data = FlatVector::GetData<uint64_t>(args.data[0]);
    auto &result_mask = FlatVector::Validity(result);

    for (idx_t i = 0; i < args.size(); i++) {
        auto block = block_data[i];

        try {
            int x, y, z;
            quadbin::cell_to_tile(block, x, y, z);

            double min_lon, min_lat, max_lon, max_lat;
            quadbin::tile_to_bbox_wgs84(x, y, z, min_lon, min_lat, max_lon, max_lat);

            auto wkb = MakeWKBPolygon(min_lon, min_lat, max_lon, max_lat);
            FlatVector::GetData<string_t>(result)[i] =
                StringVector::AddStringOrBlob(result,
                    reinterpret_cast<const char*>(wkb.data()), wkb.size());
        } catch (const std::exception &e) {
            result_mask.SetInvalid(i);
        }
    }
}

// ST_AsPolygon(block UBIGINT, band BLOB, metadata VARCHAR, band_num INTEGER) -> GEOMETRY
// Explicit band number (0-based) for selecting the band from the raster metadata.
static void STAsPolygonWithBandNumFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());
    args.data[3].Flatten(args.size());

    auto block_data = FlatVector::GetData<uint64_t>(args.data[0]);
    auto band_idx_data = FlatVector::GetData<int32_t>(args.data[3]);
    auto metadata_data = FlatVector::GetData<string_t>(args.data[2]);
    auto &result_mask = FlatVector::Validity(result);

    for (idx_t i = 0; i < args.size(); i++) {
        auto block = block_data[i];
        auto metadata_str = metadata_data[i].GetString();
        auto band_num = band_idx_data[i];

        if (band_num < 0) {
            result_mask.SetInvalid(i);
            continue;
        }

        try {
            auto meta = raquet::parse_metadata(metadata_str);

            if (band_num >= static_cast<int32_t>(meta.bands.size())) {
                throw InvalidInputException(
                    "ST_AsPolygon: band_num %d is out of range (raster has %d band(s))",
                    band_num, static_cast<int>(meta.bands.size()));
            }

            int x, y, z;
            quadbin::cell_to_tile(block, x, y, z);

            double min_lon, min_lat, max_lon, max_lat;
            quadbin::tile_to_bbox_wgs84(x, y, z, min_lon, min_lat, max_lon, max_lat);

            auto wkb = MakeWKBPolygon(min_lon, min_lat, max_lon, max_lat);
            FlatVector::GetData<string_t>(result)[i] =
                StringVector::AddStringOrBlob(result,
                    reinterpret_cast<const char*>(wkb.data()), wkb.size());
        } catch (const InvalidInputException &) {
            throw;
        } catch (const std::exception &e) {
            result_mask.SetInvalid(i);
        }
    }
}

// ============================================================================
// Function registration
// ============================================================================

void RegisterAsPolygonFunctions(ExtensionLoader &loader) {
    // ST_AsPolygon(block, metadata) -> GEOMETRY
    ScalarFunction fn("ST_AsPolygon",
        {LogicalType::UBIGINT, LogicalType::VARCHAR},
        LogicalType::GEOMETRY(),
        STAsPolygonFunction);
    loader.RegisterFunction(fn);

    // ST_AsPolygon(block, band, metadata) -> GEOMETRY
    ScalarFunction fn_band("ST_AsPolygon",
        {LogicalType::UBIGINT, LogicalType::BLOB, LogicalType::VARCHAR},
        LogicalType::GEOMETRY(),
        STAsPolygonWithBandFunction);
    loader.RegisterFunction(fn_band);

    // ST_AsPolygon(block, band, metadata, band_num) -> GEOMETRY
    ScalarFunction fn_band_num("ST_AsPolygon",
        {LogicalType::UBIGINT, LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::INTEGER},
        LogicalType::GEOMETRY(),
        STAsPolygonWithBandNumFunction);
    loader.RegisterFunction(fn_band_num);
}

} // namespace duckdb
