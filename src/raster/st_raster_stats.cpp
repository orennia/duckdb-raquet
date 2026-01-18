#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "band_decoder.hpp"
#include "quadbin.hpp"
#include <cmath>
#include <limits>

namespace duckdb {

struct RasterStats {
    int64_t count = 0;
    double sum = 0.0;
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    double mean = 0.0;
    double stddev = 0.0;

    // For Welford's online algorithm
    double m2 = 0.0;

    void add_value(double val) {
        count++;
        sum += val;
        if (val < min) min = val;
        if (val > max) max = val;

        // Welford's online algorithm for mean and variance
        double delta = val - mean;
        mean += delta / count;
        double delta2 = val - mean;
        m2 += delta * delta2;
    }

    void finalize() {
        if (count > 0) {
            mean = sum / count;
            if (count > 1) {
                stddev = std::sqrt(m2 / (count - 1));
            }
        } else {
            min = 0;
            max = 0;
        }
    }
};

// ST_RasterSummaryStats(band BLOB, dtype VARCHAR, width INT, height INT, compression VARCHAR, nodata DOUBLE)
// -> STRUCT(count BIGINT, sum DOUBLE, mean DOUBLE, min DOUBLE, max DOUBLE, stddev DOUBLE)
static void STRasterSummaryStatsFunction(DataChunk &args, ExpressionState &state, Vector &result) {
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

    auto &struct_entries = StructVector::GetEntries(result);
    auto count_data = FlatVector::GetData<int64_t>(*struct_entries[0]);
    auto sum_data = FlatVector::GetData<double>(*struct_entries[1]);
    auto mean_data = FlatVector::GetData<double>(*struct_entries[2]);
    auto min_data = FlatVector::GetData<double>(*struct_entries[3]);
    auto max_data = FlatVector::GetData<double>(*struct_entries[4]);
    auto stddev_data = FlatVector::GetData<double>(*struct_entries[5]);

    for (idx_t i = 0; i < args.size(); i++) {
        auto band = band_data[i];
        auto dtype = dtype_data[i].GetString();
        auto width = width_data[i];
        auto height = height_data[i];
        auto compression = compression_data[i].GetString();
        auto nodata = nodata_data[i];

        bool compressed = (compression == "gzip");
        bool has_nodata = !std::isnan(nodata);

        auto values = raquet::decode_band(
            reinterpret_cast<const uint8_t*>(band.GetData()),
            band.GetSize(),
            dtype, width, height, compressed
        );

        RasterStats stats;
        for (double val : values) {
            if (has_nodata && val == nodata) {
                continue;  // Skip nodata values
            }
            stats.add_value(val);
        }
        stats.finalize();

        count_data[i] = stats.count;
        sum_data[i] = stats.sum;
        mean_data[i] = stats.mean;
        min_data[i] = stats.min;
        max_data[i] = stats.max;
        stddev_data[i] = stats.stddev;
    }

    result.SetVectorType(VectorType::FLAT_VECTOR);
}

// Simplified version without nodata
static void STRasterSummaryStatsSimpleFunction(DataChunk &args, ExpressionState &state, Vector &result) {
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

    auto &struct_entries = StructVector::GetEntries(result);
    auto count_data = FlatVector::GetData<int64_t>(*struct_entries[0]);
    auto sum_data = FlatVector::GetData<double>(*struct_entries[1]);
    auto mean_data = FlatVector::GetData<double>(*struct_entries[2]);
    auto min_data = FlatVector::GetData<double>(*struct_entries[3]);
    auto max_data = FlatVector::GetData<double>(*struct_entries[4]);
    auto stddev_data = FlatVector::GetData<double>(*struct_entries[5]);

    for (idx_t i = 0; i < args.size(); i++) {
        auto band = band_data[i];
        auto dtype = dtype_data[i].GetString();
        auto width = width_data[i];
        auto height = height_data[i];
        auto compression = compression_data[i].GetString();

        bool compressed = (compression == "gzip");

        auto values = raquet::decode_band(
            reinterpret_cast<const uint8_t*>(band.GetData()),
            band.GetSize(),
            dtype, width, height, compressed
        );

        RasterStats stats;
        for (double val : values) {
            stats.add_value(val);
        }
        stats.finalize();

        count_data[i] = stats.count;
        sum_data[i] = stats.sum;
        mean_data[i] = stats.mean;
        min_data[i] = stats.min;
        max_data[i] = stats.max;
        stddev_data[i] = stats.stddev;
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
