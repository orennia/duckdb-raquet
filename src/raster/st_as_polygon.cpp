#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "band_decoder.hpp"
#include "quadbin.hpp"
#include "raquet_metadata.hpp"
#include <cmath>
#include <cstring>
#include <vector>

namespace duckdb {

// ============================================================================
// WKB geometry helpers
// ============================================================================

static void AppendUint32LE(std::vector<uint8_t> &buf, uint32_t v) {
    buf.insert(buf.end(), reinterpret_cast<uint8_t*>(&v),
               reinterpret_cast<uint8_t*>(&v) + 4);
}

static void AppendDoubleLE(std::vector<uint8_t> &buf, double v) {
    buf.insert(buf.end(), reinterpret_cast<uint8_t*>(&v),
               reinterpret_cast<uint8_t*>(&v) + 8);
}

// Build a WKB POLYGON from a WGS84 bounding box (little-endian).
// Returns a closed 5-point rectangular ring: 1 + 4 + 4 + 4 + 5*16 = 93 bytes.
static std::vector<uint8_t> MakeWKBPolygon(double min_lon, double min_lat,
                                            double max_lon, double max_lat) {
    std::vector<uint8_t> wkb;
    wkb.reserve(93);
    wkb.push_back(0x01);    // little-endian
    AppendUint32LE(wkb, 3); // Polygon
    AppendUint32LE(wkb, 1); // 1 ring
    AppendUint32LE(wkb, 5); // 5 points
    // SW, SE, NE, NW, SW (closed ring)
    AppendDoubleLE(wkb, min_lon); AppendDoubleLE(wkb, min_lat);
    AppendDoubleLE(wkb, max_lon); AppendDoubleLE(wkb, min_lat);
    AppendDoubleLE(wkb, max_lon); AppendDoubleLE(wkb, max_lat);
    AppendDoubleLE(wkb, min_lon); AppendDoubleLE(wkb, max_lat);
    AppendDoubleLE(wkb, min_lon); AppendDoubleLE(wkb, min_lat);
    return wkb;
}

// Geographic bounding box in WGS84
struct GeoRect {
    double min_lon, min_lat, max_lon, max_lat;
};

// Build a WKB MULTIPOLYGON from a list of WGS84 bounding boxes.
// Each rect becomes one rectangular polygon member.
// Layout: 9 header + N * 93 bytes per polygon.
static std::vector<uint8_t> MakeWKBMultiPolygon(const std::vector<GeoRect> &rects) {
    std::vector<uint8_t> wkb;
    wkb.reserve(9 + 93 * rects.size());
    wkb.push_back(0x01);    // little-endian
    AppendUint32LE(wkb, 6); // MultiPolygon
    AppendUint32LE(wkb, static_cast<uint32_t>(rects.size()));
    for (const auto &r : rects) {
        wkb.push_back(0x01);    // little-endian
        AppendUint32LE(wkb, 3); // Polygon
        AppendUint32LE(wkb, 1); // 1 ring
        AppendUint32LE(wkb, 5); // 5 points
        AppendDoubleLE(wkb, r.min_lon); AppendDoubleLE(wkb, r.min_lat);
        AppendDoubleLE(wkb, r.max_lon); AppendDoubleLE(wkb, r.min_lat);
        AppendDoubleLE(wkb, r.max_lon); AppendDoubleLE(wkb, r.max_lat);
        AppendDoubleLE(wkb, r.min_lon); AppendDoubleLE(wkb, r.max_lat);
        AppendDoubleLE(wkb, r.min_lon); AppendDoubleLE(wkb, r.min_lat);
    }
    return wkb;
}

// ============================================================================
// Nodata-aware polygon helpers
// ============================================================================

// Decode all pixel values for a single band, handling sequential and interleaved layouts.
static std::vector<double> DecodeBandAll(const uint8_t *data, size_t data_size,
                                          const raquet::RaquetMetadata &meta,
                                          int band_idx) {
    std::string dtype = meta.get_band_type(band_idx);
    int width = meta.block_width;
    int height = meta.block_height;
    if (meta.is_interleaved()) {
        return raquet::decode_band_interleaved(data, data_size, dtype, width, height,
                                               band_idx, meta.num_bands(),
                                               meta.compression);
    } else {
        bool compressed = (meta.compression == "gzip");
        return raquet::decode_band(data, data_size, dtype, width, height, compressed);
    }
}

// A rectangle in pixel space: [x_start, x_end) × [y_start, y_end).
struct PixelRect { int x_start, x_end, y_start, y_end; };

// Scan pixel values row by row, merging adjacent rows with identical run patterns
// into maximal rectangles of valid (non-nodata) pixels.
static std::vector<PixelRect> FindValidPixelRects(const std::vector<double> &values,
                                                   int width, int height,
                                                   const raquet::RaquetMetadata &meta,
                                                   int band_idx) {
    // Active rects being extended downward
    struct ActiveRect { int x_start, x_end, y_start; };
    std::vector<ActiveRect> active;
    std::vector<PixelRect> result;

    for (int y = 0; y < height; y++) {
        // Find contiguous runs of valid pixels in this row
        std::vector<std::pair<int,int>> runs;
        int run_start = -1;
        for (int x = 0; x < width; x++) {
            bool valid = !meta.is_nodata(band_idx, values[y * width + x]);
            if (valid) {
                if (run_start < 0) run_start = x;
            } else {
                if (run_start >= 0) { runs.push_back({run_start, x}); run_start = -1; }
            }
        }
        if (run_start >= 0) runs.push_back({run_start, width});

        // Match each run to an existing active rect with the same x range.
        // Use a 'used' flag to avoid O(n) vector erase.
        std::vector<bool> used(active.size(), false);
        std::vector<ActiveRect> new_active;
        for (const auto &run : runs) {
            bool matched = false;
            for (size_t j = 0; j < active.size(); j++) {
                if (!used[j] && active[j].x_start == run.first && active[j].x_end == run.second) {
                    new_active.push_back(active[j]);
                    used[j] = true;
                    matched = true;
                    break;
                }
            }
            if (!matched) new_active.push_back({run.first, run.second, y});
        }

        // Any active rects not matched are completed at this row
        for (size_t j = 0; j < active.size(); j++) {
            if (!used[j]) result.push_back({active[j].x_start, active[j].x_end, active[j].y_start, y});
        }
        active = std::move(new_active);
    }

    // Flush remaining active rects
    for (const auto &r : active) result.push_back({r.x_start, r.x_end, r.y_start, height});
    return result;
}

// Convert a pixel-space rectangle to a WGS84 geographic bounding box.
// Pixels are linearly spaced in Mercator; each pixel maps to a Mercator rectangle
// whose corners are then individually projected to WGS84 (lon/lat).
//
// Coordinate system:
//   - pixel row 0 is at the TOP of the tile (highest Mercator y = max_y)
//   - row increases downward (decreasing Mercator y)
//   - pixel col 0 is at the LEFT of the tile (lowest Mercator x = min_x)
static GeoRect PixelRectToGeo(const PixelRect &rect, int block_width, int block_height,
                               double merc_min_x, double merc_min_y,
                               double merc_max_x, double merc_max_y) {
    double pixel_w = (merc_max_x - merc_min_x) / block_width;
    double pixel_h = (merc_max_y - merc_min_y) / block_height;

    double mx0 = merc_min_x + rect.x_start * pixel_w;
    double mx1 = merc_min_x + rect.x_end   * pixel_w;
    double my0 = merc_max_y - rect.y_end   * pixel_h; // south edge of rect
    double my1 = merc_max_y - rect.y_start * pixel_h; // north edge of rect

    // Mercator → WGS84: lon = x/R * (180/π), lat = 2*atan(exp(y/R)) * (180/π) - 90
    const double R = quadbin::EARTH_RADIUS;
    const double DEG = 180.0 / quadbin::PI;
    GeoRect g;
    g.min_lon = mx0 / R * DEG;
    g.max_lon = mx1 / R * DEG;
    g.min_lat = 2.0 * std::atan(std::exp(my0 / R)) * DEG - 90.0;
    g.max_lat = 2.0 * std::atan(std::exp(my1 / R)) * DEG - 90.0;
    return g;
}

// Core of ST_AsPolygon when band data is available.
// Returns the WKB geometry (POLYGON for full-tile case, MULTIPOLYGON for filtered case),
// or an empty vector if all pixels are nodata.
static std::vector<uint8_t> ComputeBandPolygon(uint64_t block, const string_t &band,
                                                const raquet::RaquetMetadata &meta,
                                                int band_idx) {
    // Fast path: no nodata defined for this band → return full tile bounding box as POLYGON
    bool has_nodata = (band_idx < static_cast<int>(meta.band_info.size()) &&
                       meta.band_info[band_idx].has_nodata);
    if (!has_nodata) {
        int tx, ty, tz;
        quadbin::cell_to_tile(block, tx, ty, tz);
        double min_lon, min_lat, max_lon, max_lat;
        quadbin::tile_to_bbox_wgs84(tx, ty, tz, min_lon, min_lat, max_lon, max_lat);
        return MakeWKBPolygon(min_lon, min_lat, max_lon, max_lat);
    }

    // Decode all pixels for this band
    auto values = DecodeBandAll(reinterpret_cast<const uint8_t*>(band.GetData()),
                                band.GetSize(), meta, band_idx);

    // Find valid pixel rectangles (excludes nodata pixels)
    auto pixel_rects = FindValidPixelRects(values, meta.block_width, meta.block_height,
                                           meta, band_idx);
    if (pixel_rects.empty()) return {}; // all pixels are nodata → NULL

    // Get tile bounds in Mercator for coordinate conversion
    int tx, ty, tz;
    quadbin::cell_to_tile(block, tx, ty, tz);
    double merc_min_x, merc_min_y, merc_max_x, merc_max_y;
    quadbin::tile_to_bbox_mercator(tx, ty, tz, merc_min_x, merc_min_y, merc_max_x, merc_max_y);

    // Convert pixel rectangles to WGS84 geographic bounds
    std::vector<GeoRect> geo_rects;
    geo_rects.reserve(pixel_rects.size());
    for (const auto &pr : pixel_rects) {
        geo_rects.push_back(PixelRectToGeo(pr, meta.block_width, meta.block_height,
                                           merc_min_x, merc_min_y, merc_max_x, merc_max_y));
    }

    return MakeWKBMultiPolygon(geo_rects);
}

// ============================================================================
// ST_AsPolygon scalar functions
// ============================================================================

// ST_AsPolygon(block UBIGINT, metadata VARCHAR) -> GEOMETRY
// No band data available; returns the full tile bounding polygon.
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
        } catch (const std::exception &) {
            result_mask.SetInvalid(i);
        }
    }
}

// ST_AsPolygon(block UBIGINT, band BLOB, metadata VARCHAR) -> GEOMETRY
// Uses band 0 (first band) to exclude nodata pixels from the output polygon.
// Returns NULL if all pixels are nodata.
static void STAsPolygonWithBandFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());

    auto block_data    = FlatVector::GetData<uint64_t>(args.data[0]);
    auto band_data     = FlatVector::GetData<string_t>(args.data[1]);
    auto metadata_data = FlatVector::GetData<string_t>(args.data[2]);
    auto &result_mask  = FlatVector::Validity(result);

    for (idx_t i = 0; i < args.size(); i++) {
        auto block        = block_data[i];
        auto band         = band_data[i];
        auto metadata_str = metadata_data[i].GetString();

        if (band.GetSize() == 0) {
            result_mask.SetInvalid(i);
            continue;
        }

        try {
            auto meta = raquet::parse_metadata(metadata_str);
            auto wkb  = ComputeBandPolygon(block, band, meta, 0);
            if (wkb.empty()) {
                result_mask.SetInvalid(i); // all pixels nodata
            } else {
                FlatVector::GetData<string_t>(result)[i] =
                    StringVector::AddStringOrBlob(result,
                        reinterpret_cast<const char*>(wkb.data()), wkb.size());
            }
        } catch (const std::exception &) {
            result_mask.SetInvalid(i);
        }
    }
}

// ST_AsPolygon(block UBIGINT, band BLOB, metadata VARCHAR, band_num INTEGER) -> GEOMETRY
// Uses the specified band (0-based) to exclude nodata pixels from the output polygon.
// Returns NULL if all pixels are nodata or band_num is negative.
static void STAsPolygonWithBandNumFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());
    args.data[3].Flatten(args.size());

    auto block_data    = FlatVector::GetData<uint64_t>(args.data[0]);
    auto band_data     = FlatVector::GetData<string_t>(args.data[1]);
    auto metadata_data = FlatVector::GetData<string_t>(args.data[2]);
    auto band_idx_data = FlatVector::GetData<int32_t>(args.data[3]);
    auto &result_mask  = FlatVector::Validity(result);

    for (idx_t i = 0; i < args.size(); i++) {
        auto block        = block_data[i];
        auto band         = band_data[i];
        auto metadata_str = metadata_data[i].GetString();
        auto band_num     = band_idx_data[i];

        if (band_num < 0) {
            result_mask.SetInvalid(i);
            continue;
        }

        if (band.GetSize() == 0) {
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

            auto wkb = ComputeBandPolygon(block, band, meta, band_num);
            if (wkb.empty()) {
                result_mask.SetInvalid(i); // all pixels nodata
            } else {
                FlatVector::GetData<string_t>(result)[i] =
                    StringVector::AddStringOrBlob(result,
                        reinterpret_cast<const char*>(wkb.data()), wkb.size());
            }
        } catch (const InvalidInputException &) {
            throw;
        } catch (const std::exception &) {
            result_mask.SetInvalid(i);
        }
    }
}

// ============================================================================
// Function registration
// ============================================================================

void RegisterAsPolygonFunctions(ExtensionLoader &loader) {
    // ST_AsPolygon(block, metadata) -> GEOMETRY  (no band data; returns full tile bbox)
    ScalarFunction fn("ST_AsPolygon",
        {LogicalType::UBIGINT, LogicalType::VARCHAR},
        LogicalType::GEOMETRY(),
        STAsPolygonFunction);
    loader.RegisterFunction(fn);

    // ST_AsPolygon(block, band, metadata) -> GEOMETRY  (uses band 0 for nodata filtering)
    ScalarFunction fn_band("ST_AsPolygon",
        {LogicalType::UBIGINT, LogicalType::BLOB, LogicalType::VARCHAR},
        LogicalType::GEOMETRY(),
        STAsPolygonWithBandFunction);
    loader.RegisterFunction(fn_band);

    // ST_AsPolygon(block, band, metadata, band_num) -> GEOMETRY  (explicit 0-based band index)
    ScalarFunction fn_band_num("ST_AsPolygon",
        {LogicalType::UBIGINT, LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::INTEGER},
        LogicalType::GEOMETRY(),
        STAsPolygonWithBandNumFunction);
    loader.RegisterFunction(fn_band_num);
}

} // namespace duckdb
