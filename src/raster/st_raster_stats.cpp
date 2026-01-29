#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "band_decoder.hpp"
#include "quadbin.hpp"
#include <cmath>

namespace duckdb {

// ST_RasterSummaryStats(band BLOB, dtype VARCHAR, width INT, height INT, compression VARCHAR, nodata DOUBLE)
// -> STRUCT(count BIGINT, sum DOUBLE, mean DOUBLE, min DOUBLE, max DOUBLE, stddev DOUBLE)
static void STRasterSummaryStatsFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    // Flatten all input vectors to ensure consistent access
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());
    args.data[3].Flatten(args.size());
    args.data[4].Flatten(args.size());
    args.data[5].Flatten(args.size());

    auto &band_vec = args.data[0];
    auto &dtype_vec = args.data[1];
    auto &width_vec = args.data[2];
    auto &height_vec = args.data[3];
    auto &compression_vec = args.data[4];
    auto &nodata_vec = args.data[5];

    auto band_data = FlatVector::GetData<string_t>(band_vec);
    auto dtype_data = FlatVector::GetData<string_t>(dtype_vec);
    auto width_data = FlatVector::GetData<int32_t>(width_vec);
    auto height_data = FlatVector::GetData<int32_t>(height_vec);
    auto compression_data = FlatVector::GetData<string_t>(compression_vec);
    auto nodata_data = FlatVector::GetData<double>(nodata_vec);

    auto &band_validity = FlatVector::Validity(band_vec);
    auto &nodata_validity = FlatVector::Validity(nodata_vec);

    auto &struct_entries = StructVector::GetEntries(result);
    auto count_data = FlatVector::GetData<int64_t>(*struct_entries[0]);
    auto sum_data = FlatVector::GetData<double>(*struct_entries[1]);
    auto mean_data = FlatVector::GetData<double>(*struct_entries[2]);
    auto min_data = FlatVector::GetData<double>(*struct_entries[3]);
    auto max_data = FlatVector::GetData<double>(*struct_entries[4]);
    auto stddev_data = FlatVector::GetData<double>(*struct_entries[5]);

    auto &result_validity = FlatVector::Validity(result);

    for (idx_t i = 0; i < args.size(); i++) {
        // Check for NULL band data
        if (!band_validity.RowIsValid(i)) {
            result_validity.SetInvalid(i);
            continue;
        }

        auto band = band_data[i];
        if (band.GetSize() == 0) {
            result_validity.SetInvalid(i);
            continue;
        }

        auto dtype = dtype_data[i].GetString();
        auto width = width_data[i];
        auto height = height_data[i];
        auto compression = compression_data[i].GetString();
        auto nodata = nodata_data[i];

        bool compressed = (compression == "gzip");
        // Support NaN as a valid nodata value (Zarr v3 convention)
        bool has_nodata = nodata_validity.RowIsValid(i);

        try {
            // Use streaming stats for better performance (avoids allocating full pixel array)
            auto stats = raquet::compute_band_stats(
                reinterpret_cast<const uint8_t*>(band.GetData()),
                band.GetSize(),
                dtype, width, height, compressed,
                has_nodata, nodata
            );

            count_data[i] = stats.count;
            sum_data[i] = stats.sum;
            mean_data[i] = stats.mean;
            min_data[i] = stats.min;
            max_data[i] = stats.max;
            stddev_data[i] = stats.stddev;
        } catch (const std::exception &e) {
            // On decompression failure, return NULL
            result_validity.SetInvalid(i);
        }
    }

    result.SetVectorType(VectorType::FLAT_VECTOR);
}

// Simplified version without nodata
static void STRasterSummaryStatsSimpleFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    // Flatten all input vectors to ensure consistent access
    args.data[0].Flatten(args.size());
    args.data[1].Flatten(args.size());
    args.data[2].Flatten(args.size());
    args.data[3].Flatten(args.size());
    args.data[4].Flatten(args.size());

    auto &band_vec = args.data[0];
    auto &dtype_vec = args.data[1];
    auto &width_vec = args.data[2];
    auto &height_vec = args.data[3];
    auto &compression_vec = args.data[4];

    auto band_data = FlatVector::GetData<string_t>(band_vec);
    auto dtype_data = FlatVector::GetData<string_t>(dtype_vec);
    auto width_data = FlatVector::GetData<int32_t>(width_vec);
    auto height_data = FlatVector::GetData<int32_t>(height_vec);
    auto compression_data = FlatVector::GetData<string_t>(compression_vec);

    auto &band_validity = FlatVector::Validity(band_vec);

    auto &struct_entries = StructVector::GetEntries(result);
    auto count_data = FlatVector::GetData<int64_t>(*struct_entries[0]);
    auto sum_data = FlatVector::GetData<double>(*struct_entries[1]);
    auto mean_data = FlatVector::GetData<double>(*struct_entries[2]);
    auto min_data = FlatVector::GetData<double>(*struct_entries[3]);
    auto max_data = FlatVector::GetData<double>(*struct_entries[4]);
    auto stddev_data = FlatVector::GetData<double>(*struct_entries[5]);

    auto &result_validity = FlatVector::Validity(result);

    for (idx_t i = 0; i < args.size(); i++) {
        // Check for NULL band data
        if (!band_validity.RowIsValid(i)) {
            result_validity.SetInvalid(i);
            continue;
        }

        auto band = band_data[i];
        if (band.GetSize() == 0) {
            result_validity.SetInvalid(i);
            continue;
        }

        auto dtype = dtype_data[i].GetString();
        auto width = width_data[i];
        auto height = height_data[i];
        auto compression = compression_data[i].GetString();

        bool compressed = (compression == "gzip");

        try {
            // Use streaming stats for better performance (avoids allocating full pixel array)
            auto stats = raquet::compute_band_stats(
                reinterpret_cast<const uint8_t*>(band.GetData()),
                band.GetSize(),
                dtype, width, height, compressed,
                false, 0.0  // no nodata filtering
            );

            count_data[i] = stats.count;
            sum_data[i] = stats.sum;
            mean_data[i] = stats.mean;
            min_data[i] = stats.min;
            max_data[i] = stats.max;
            stddev_data[i] = stats.stddev;
        } catch (const std::exception &e) {
            // On decompression failure, return NULL
            result_validity.SetInvalid(i);
        }
    }

    result.SetVectorType(VectorType::FLAT_VECTOR);
}

void RegisterRasterStatsFunctions(ExtensionLoader &loader) {
    // Define the stats struct type
    child_list_t<LogicalType> stats_struct;
    stats_struct.push_back(make_pair("count", LogicalType::BIGINT));
    stats_struct.push_back(make_pair("sum", LogicalType::DOUBLE));
    stats_struct.push_back(make_pair("mean", LogicalType::DOUBLE));
    stats_struct.push_back(make_pair("min", LogicalType::DOUBLE));
    stats_struct.push_back(make_pair("max", LogicalType::DOUBLE));
    stats_struct.push_back(make_pair("stddev", LogicalType::DOUBLE));
    auto stats_type = LogicalType::STRUCT(stats_struct);

    // ST_RasterSummaryStats with nodata
    ScalarFunction stats_fn("ST_RasterSummaryStats",
        {LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::INTEGER,
         LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::DOUBLE},
        stats_type,
        STRasterSummaryStatsFunction);
    loader.RegisterFunction(stats_fn);

    // ST_RasterSummaryStats without nodata
    ScalarFunction stats_simple_fn("ST_RasterSummaryStats",
        {LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::INTEGER,
         LogicalType::INTEGER, LogicalType::VARCHAR},
        stats_type,
        STRasterSummaryStatsSimpleFunction);
    loader.RegisterFunction(stats_simple_fn);
}

} // namespace duckdb
