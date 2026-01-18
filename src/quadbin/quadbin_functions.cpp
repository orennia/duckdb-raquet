#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "quadbin.hpp"

namespace duckdb {

// quadbin_from_tile(x, y, z) -> UBIGINT
static void QuadbinFromTileFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &x_vec = args.data[0];
    auto &y_vec = args.data[1];
    auto &z_vec = args.data[2];

    TernaryExecutor::Execute<int32_t, int32_t, int32_t, uint64_t>(
        x_vec, y_vec, z_vec, result, args.size(),
        [&](int32_t x, int32_t y, int32_t z) {
            return quadbin::tile_to_cell(x, y, z);
        });
}

// quadbin_to_tile(cell) -> STRUCT(x INT, y INT, z INT)
static void QuadbinToTileFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &cell_vec = args.data[0];

    // Result is a struct with x, y, z fields
    auto &struct_entries = StructVector::GetEntries(result);
    auto &x_result = *struct_entries[0];
    auto &y_result = *struct_entries[1];
    auto &z_result = *struct_entries[2];

    UnaryExecutor::Execute<uint64_t, int32_t>(cell_vec, x_result, args.size(),
        [](uint64_t cell) {
            int x, y, z;
            quadbin::cell_to_tile(cell, x, y, z);
            return x;
        });

    UnaryExecutor::Execute<uint64_t, int32_t>(cell_vec, y_result, args.size(),
        [](uint64_t cell) {
            int x, y, z;
            quadbin::cell_to_tile(cell, x, y, z);
            return y;
        });

    UnaryExecutor::Execute<uint64_t, int32_t>(cell_vec, z_result, args.size(),
        [](uint64_t cell) {
            int x, y, z;
            quadbin::cell_to_tile(cell, x, y, z);
            return z;
        });

    result.SetVectorType(VectorType::FLAT_VECTOR);
}

// quadbin_from_lonlat(lon, lat, resolution) -> UBIGINT
static void QuadbinFromLonLatFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &lon_vec = args.data[0];
    auto &lat_vec = args.data[1];
    auto &res_vec = args.data[2];

    TernaryExecutor::Execute<double, double, int32_t, uint64_t>(
        lon_vec, lat_vec, res_vec, result, args.size(),
        [](double lon, double lat, int32_t resolution) {
            return quadbin::lonlat_to_cell(lon, lat, resolution);
        });
}

// quadbin_to_lonlat(cell) -> STRUCT(lon DOUBLE, lat DOUBLE)
static void QuadbinToLonLatFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &cell_vec = args.data[0];

    auto &struct_entries = StructVector::GetEntries(result);
    auto &lon_result = *struct_entries[0];
    auto &lat_result = *struct_entries[1];

    UnaryExecutor::Execute<uint64_t, double>(cell_vec, lon_result, args.size(),
        [](uint64_t cell) {
            double lon, lat;
            quadbin::cell_to_lonlat(cell, lon, lat);
            return lon;
        });

    UnaryExecutor::Execute<uint64_t, double>(cell_vec, lat_result, args.size(),
        [](uint64_t cell) {
            double lon, lat;
            quadbin::cell_to_lonlat(cell, lon, lat);
            return lat;
        });

    result.SetVectorType(VectorType::FLAT_VECTOR);
}

// quadbin_resolution(cell) -> INT
static void QuadbinResolutionFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &cell_vec = args.data[0];

    UnaryExecutor::Execute<uint64_t, int32_t>(cell_vec, result, args.size(),
        [](uint64_t cell) {
            return quadbin::cell_to_resolution(cell);
        });
}

// quadbin_to_bbox(cell) -> STRUCT(min_lon DOUBLE, min_lat DOUBLE, max_lon DOUBLE, max_lat DOUBLE)
static void QuadbinToBboxFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &cell_vec = args.data[0];

    auto &struct_entries = StructVector::GetEntries(result);
    auto &min_lon_result = *struct_entries[0];
    auto &min_lat_result = *struct_entries[1];
    auto &max_lon_result = *struct_entries[2];
    auto &max_lat_result = *struct_entries[3];

    for (idx_t i = 0; i < args.size(); i++) {
        auto cell = FlatVector::GetData<uint64_t>(cell_vec)[i];

        int x, y, z;
        quadbin::cell_to_tile(cell, x, y, z);

        double min_lon, min_lat, max_lon, max_lat;
        quadbin::tile_to_bbox_wgs84(x, y, z, min_lon, min_lat, max_lon, max_lat);

        FlatVector::GetData<double>(min_lon_result)[i] = min_lon;
        FlatVector::GetData<double>(min_lat_result)[i] = min_lat;
        FlatVector::GetData<double>(max_lon_result)[i] = max_lon;
        FlatVector::GetData<double>(max_lat_result)[i] = max_lat;
    }

    result.SetVectorType(VectorType::FLAT_VECTOR);
}

// quadbin_pixel_xy(lon, lat, resolution, tile_size) -> STRUCT(pixel_x INT, pixel_y INT)
static void QuadbinPixelXYFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &lon_vec = args.data[0];
    auto &lat_vec = args.data[1];
    auto &res_vec = args.data[2];
    auto &size_vec = args.data[3];

    auto &struct_entries = StructVector::GetEntries(result);
    auto &px_result = *struct_entries[0];
    auto &py_result = *struct_entries[1];

    for (idx_t i = 0; i < args.size(); i++) {
        auto lon = FlatVector::GetData<double>(lon_vec)[i];
        auto lat = FlatVector::GetData<double>(lat_vec)[i];
        auto resolution = FlatVector::GetData<int32_t>(res_vec)[i];
        auto tile_size = FlatVector::GetData<int32_t>(size_vec)[i];

        int pixel_x, pixel_y, tile_x, tile_y;
        quadbin::lonlat_to_pixel(lon, lat, resolution, tile_size, pixel_x, pixel_y, tile_x, tile_y);

        FlatVector::GetData<int32_t>(px_result)[i] = pixel_x;
        FlatVector::GetData<int32_t>(py_result)[i] = pixel_y;
    }

    result.SetVectorType(VectorType::FLAT_VECTOR);
}

void RegisterQuadbinFunctions(ExtensionLoader &loader) {
    // quadbin_from_tile(x, y, z) -> UBIGINT
    ScalarFunction from_tile("quadbin_from_tile",
        {LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER},
        LogicalType::UBIGINT,
        QuadbinFromTileFunction);
    loader.RegisterFunction(from_tile);

    // quadbin_to_tile(cell) -> STRUCT
    child_list_t<LogicalType> tile_struct;
    tile_struct.push_back(make_pair("x", LogicalType::INTEGER));
    tile_struct.push_back(make_pair("y", LogicalType::INTEGER));
    tile_struct.push_back(make_pair("z", LogicalType::INTEGER));

    ScalarFunction to_tile("quadbin_to_tile",
        {LogicalType::UBIGINT},
        LogicalType::STRUCT(tile_struct),
        QuadbinToTileFunction);
    loader.RegisterFunction(to_tile);

    // quadbin_from_lonlat(lon, lat, resolution) -> UBIGINT
    ScalarFunction from_lonlat("quadbin_from_lonlat",
        {LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::INTEGER},
        LogicalType::UBIGINT,
        QuadbinFromLonLatFunction);
    loader.RegisterFunction(from_lonlat);

    // quadbin_to_lonlat(cell) -> STRUCT
    child_list_t<LogicalType> lonlat_struct;
    lonlat_struct.push_back(make_pair("lon", LogicalType::DOUBLE));
    lonlat_struct.push_back(make_pair("lat", LogicalType::DOUBLE));

    ScalarFunction to_lonlat("quadbin_to_lonlat",
        {LogicalType::UBIGINT},
        LogicalType::STRUCT(lonlat_struct),
        QuadbinToLonLatFunction);
    loader.RegisterFunction(to_lonlat);

    // quadbin_resolution(cell) -> INT
    ScalarFunction resolution("quadbin_resolution",
        {LogicalType::UBIGINT},
        LogicalType::INTEGER,
        QuadbinResolutionFunction);
    loader.RegisterFunction(resolution);

    // quadbin_to_bbox(cell) -> STRUCT(min_lon, min_lat, max_lon, max_lat)
    child_list_t<LogicalType> bbox_struct;
    bbox_struct.push_back(make_pair("min_lon", LogicalType::DOUBLE));
    bbox_struct.push_back(make_pair("min_lat", LogicalType::DOUBLE));
    bbox_struct.push_back(make_pair("max_lon", LogicalType::DOUBLE));
    bbox_struct.push_back(make_pair("max_lat", LogicalType::DOUBLE));

    ScalarFunction to_bbox("quadbin_to_bbox",
        {LogicalType::UBIGINT},
        LogicalType::STRUCT(bbox_struct),
        QuadbinToBboxFunction);
    loader.RegisterFunction(to_bbox);

    // quadbin_pixel_xy(lon, lat, resolution, tile_size) -> STRUCT(pixel_x, pixel_y)
    child_list_t<LogicalType> pixel_struct;
    pixel_struct.push_back(make_pair("pixel_x", LogicalType::INTEGER));
    pixel_struct.push_back(make_pair("pixel_y", LogicalType::INTEGER));

    ScalarFunction pixel_xy("quadbin_pixel_xy",
        {LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::INTEGER, LogicalType::INTEGER},
        LogicalType::STRUCT(pixel_struct),
        QuadbinPixelXYFunction);
    loader.RegisterFunction(pixel_xy);
}

} // namespace duckdb
