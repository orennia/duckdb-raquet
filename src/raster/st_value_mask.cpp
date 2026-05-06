#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "band_decoder.hpp"
#include "band_encoder.hpp"
#include "yyjson.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace duckdb_yyjson;

namespace duckdb {

namespace {

struct MaskTileMetadata {
	std::string compression {"gzip"};
	std::string band_layout {"sequential"};
	int block_width {256};
	int block_height {256};
	std::string dtype_str {"uint8"};

	static MaskTileMetadata Parse(const std::string &json_str) {
		MaskTileMetadata meta;
		yyjson_doc *doc = yyjson_read(json_str.c_str(), json_str.length(), 0);
		if (!doc) {
			throw InvalidInputException("ST_MaskBandValues*: Invalid metadata JSON");
		}
		yyjson_val *root = yyjson_doc_get_root(doc);
		if (!yyjson_is_obj(root)) {
			yyjson_doc_free(doc);
			throw InvalidInputException("ST_MaskBandValues*: Metadata must be a JSON object");
		}

		yyjson_val *compression_val = yyjson_obj_get(root, "compression");
		if (compression_val && yyjson_is_str(compression_val)) {
			meta.compression = yyjson_get_str(compression_val);
		}
		yyjson_val *layout_val = yyjson_obj_get(root, "band_layout");
		if (layout_val && yyjson_is_str(layout_val)) {
			meta.band_layout = yyjson_get_str(layout_val);
		}

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

		yyjson_val *bands_val = yyjson_obj_get(root, "bands");
		if (bands_val && yyjson_is_arr(bands_val)) {
			size_t idx, max;
			yyjson_val *band;
			yyjson_arr_foreach(bands_val, idx, max, band) {
				if (yyjson_is_obj(band)) {
					yyjson_val *dtype_val = yyjson_obj_get(band, "type");
					if (dtype_val && yyjson_is_str(dtype_val)) {
						meta.dtype_str = yyjson_get_str(dtype_val);
						break;
					}
				}
			}
		}

		yyjson_doc_free(doc);
		return meta;
	}

	void RequireGzipSequentialBand() const {
		if (band_layout != "sequential") {
			throw InvalidInputException(
			    "ST_MaskBandValues*: only band_layout='sequential' per-band tiles are supported (got '%s')",
			    band_layout.c_str());
		}
		if (compression != "gzip" && compression != "none" && !compression.empty()) {
			throw InvalidInputException(
			    "ST_MaskBandValues*: only compression 'gzip' or 'none' is supported for sequential bands (got '%s')",
			    compression.c_str());
		}
	}
};

static bool double_in_sorted_unique(const std::vector<double> &sorted, double v) {
	return std::binary_search(sorted.begin(), sorted.end(), v);
}

static std::vector<double> DecodePixels(const string_t &band, const MaskTileMetadata &meta) {
	const bool compressed = (meta.compression == "gzip");
	return raquet::decode_band(reinterpret_cast<const uint8_t *>(band.GetData()), band.GetSize(), meta.dtype_str,
	                           meta.block_width, meta.block_height, compressed);
}

static string_t EncodeBlob(Vector &result, const std::vector<uint8_t> &raw, const MaskTileMetadata &meta) {
	std::vector<uint8_t> out;
	if (meta.compression == "gzip") {
		out = raquet::compress_gzip(raw.data(), raw.size());
	} else {
		out = raw;
	}
	return StringVector::AddStringOrBlob(result, reinterpret_cast<const char *>(out.data()), out.size());
}

static void MaskBandValuesListFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	args.data[0].Flatten(args.size());
	args.data[1].Flatten(args.size());
	args.data[2].Flatten(args.size());
	args.data[3].Flatten(args.size());

	auto band_data = FlatVector::GetData<string_t>(args.data[0]);
	auto metadata_data = FlatVector::GetData<string_t>(args.data[1]);
	auto nodata_data = FlatVector::GetData<double>(args.data[2]);
	auto list_entries = ListVector::GetData(args.data[3]);
	auto &list_child = ListVector::GetEntry(args.data[3]);
	auto child_data = FlatVector::GetData<double>(list_child);

	auto &band_validity = FlatVector::Validity(args.data[0]);
	auto &meta_validity = FlatVector::Validity(args.data[1]);
	auto &list_validity = FlatVector::Validity(args.data[3]);

	auto result_data = FlatVector::GetData<string_t>(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!band_validity.RowIsValid(i) || !meta_validity.RowIsValid(i) || !list_validity.RowIsValid(i)) {
			FlatVector::Validity(result).SetInvalid(i);
			continue;
		}
		const auto &le = list_entries[i];
		if (le.length == 0) {
			throw InvalidInputException(
			    "ST_MaskBandValuesList: allowed-values list must be non-empty (inclusion / whitelist)");
		}
		std::vector<double> allowed;
		allowed.reserve(le.length);
		for (idx_t j = 0; j < le.length; j++) {
			allowed.push_back(child_data[le.offset + j]);
		}
		std::sort(allowed.begin(), allowed.end());
		allowed.erase(std::unique(allowed.begin(), allowed.end()), allowed.end());

		const auto band = band_data[i];
		if (band.GetSize() == 0) {
			FlatVector::Validity(result).SetInvalid(i);
			continue;
		}
		auto meta = MaskTileMetadata::Parse(metadata_data[i].GetString());
		meta.RequireGzipSequentialBand();
		auto dtype = raquet::parse_dtype(meta.dtype_str);
		auto pixels = DecodePixels(band, meta);
		const double nodata = nodata_data[i];
		// Inclusion: keep pixels whose value is in the list; nodata outside the set (and non-finite).
		for (double &p : pixels) {
			if (!std::isfinite(p) || !double_in_sorted_unique(allowed, p)) {
				p = nodata;
			}
		}
		auto raw = raquet::encode_band_from_doubles(pixels, dtype);
		result_data[i] = EncodeBlob(result, raw, meta);
	}
}

static void MaskBandValuesRangeFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	args.data[0].Flatten(args.size());
	args.data[1].Flatten(args.size());
	args.data[2].Flatten(args.size());
	args.data[3].Flatten(args.size());
	args.data[4].Flatten(args.size());

	auto band_data = FlatVector::GetData<string_t>(args.data[0]);
	auto metadata_data = FlatVector::GetData<string_t>(args.data[1]);
	auto nodata_data = FlatVector::GetData<double>(args.data[2]);
	auto low_data = FlatVector::GetData<double>(args.data[3]);
	auto high_data = FlatVector::GetData<double>(args.data[4]);

	auto &band_validity = FlatVector::Validity(args.data[0]);
	auto &meta_validity = FlatVector::Validity(args.data[1]);
	auto &low_validity = FlatVector::Validity(args.data[3]);
	auto &high_validity = FlatVector::Validity(args.data[4]);

	auto result_data = FlatVector::GetData<string_t>(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!band_validity.RowIsValid(i) || !meta_validity.RowIsValid(i)) {
			FlatVector::Validity(result).SetInvalid(i);
			continue;
		}
		const auto band = band_data[i];
		if (band.GetSize() == 0) {
			FlatVector::Validity(result).SetInvalid(i);
			continue;
		}
		const bool low_set = low_validity.RowIsValid(i);
		const bool high_set = high_validity.RowIsValid(i);
		if (!low_set && !high_set) {
			throw InvalidInputException(
			    "ST_MaskBandValuesRange: at least one of low or high must be set (non-NULL)");
		}
		const double low = low_data[i];
		const double high = high_data[i];
		if (low_set && high_set && !(low <= high)) {
			throw InvalidInputException("ST_MaskBandValuesRange: low must be <= high when both are set");
		}
		auto meta = MaskTileMetadata::Parse(metadata_data[i].GetString());
		meta.RequireGzipSequentialBand();
		auto dtype = raquet::parse_dtype(meta.dtype_str);
		auto pixels = DecodePixels(band, meta);
		const double nodata = nodata_data[i];
		// Keep band: inclusive on any bound that is set. NULL low → no lower cap; NULL high → no upper cap.
		// Examples: keep slope ≥ 5, no max — ST_MaskBandValuesRange(band, metadata, nodata, 5, NULL).
		//           keep slope ≤ 5 only — ST_MaskBandValuesRange(band, metadata, nodata, NULL, 5).
		for (double &p : pixels) {
			if (!std::isfinite(p)) {
				p = nodata;
				continue;
			}
			if (low_set && p < low) {
				p = nodata;
			} else if (high_set && p > high) {
				p = nodata;
			}
		}
		auto raw = raquet::encode_band_from_doubles(pixels, dtype);
		result_data[i] = EncodeBlob(result, raw, meta);
	}
}

static void MaskBandValuesEqFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	args.data[0].Flatten(args.size());
	args.data[1].Flatten(args.size());
	args.data[2].Flatten(args.size());
	args.data[3].Flatten(args.size());

	auto band_data = FlatVector::GetData<string_t>(args.data[0]);
	auto metadata_data = FlatVector::GetData<string_t>(args.data[1]);
	auto nodata_data = FlatVector::GetData<double>(args.data[2]);
	auto value_data = FlatVector::GetData<double>(args.data[3]);

	auto &band_validity = FlatVector::Validity(args.data[0]);
	auto &meta_validity = FlatVector::Validity(args.data[1]);

	auto result_data = FlatVector::GetData<string_t>(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!band_validity.RowIsValid(i) || !meta_validity.RowIsValid(i)) {
			FlatVector::Validity(result).SetInvalid(i);
			continue;
		}
		const auto band = band_data[i];
		if (band.GetSize() == 0) {
			FlatVector::Validity(result).SetInvalid(i);
			continue;
		}
		const double eq = value_data[i];
		auto meta = MaskTileMetadata::Parse(metadata_data[i].GetString());
		meta.RequireGzipSequentialBand();
		auto dtype = raquet::parse_dtype(meta.dtype_str);
		auto pixels = DecodePixels(band, meta);
		const double nodata = nodata_data[i];
		// Inclusion: keep pixels equal to `value` only; nodata elsewhere (including non-finite).
		for (double &p : pixels) {
			if (!std::isfinite(p) || p != eq) {
				p = nodata;
			}
		}
		auto raw = raquet::encode_band_from_doubles(pixels, dtype);
		result_data[i] = EncodeBlob(result, raw, meta);
	}
}

} // namespace

void RegisterValueMaskFunctions(ExtensionLoader &loader) {
	// NULL in / NULL out would make NULL low or NULL high yield a NULL blob without running the callback;
	// SPECIAL_HANDLING lets optional bounds mean "open" on that side.
	// ST_MaskBandValuesList: inclusion — keep pixels whose value appears in the list; nodata otherwise.
	ScalarFunction list_fn(
	    "ST_MaskBandValuesList",
	    {LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::DOUBLE, LogicalType::LIST(LogicalType::DOUBLE)},
	    LogicalType::BLOB, MaskBandValuesListFunction, nullptr, nullptr, nullptr, nullptr,
	    LogicalType(LogicalTypeId::INVALID), FunctionStability::CONSISTENT, FunctionNullHandling::SPECIAL_HANDLING);
	loader.RegisterFunction(list_fn);

	// ST_MaskBandValuesRange: inclusion — keep inclusive within any bound provided; NULL low or NULL high = open on that side.
	ScalarFunction range_fn("ST_MaskBandValuesRange",
	                        {LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::DOUBLE, LogicalType::DOUBLE,
	                         LogicalType::DOUBLE},
	                        LogicalType::BLOB, MaskBandValuesRangeFunction, nullptr, nullptr, nullptr, nullptr,
	                        LogicalType(LogicalTypeId::INVALID), FunctionStability::CONSISTENT,
	                        FunctionNullHandling::SPECIAL_HANDLING);
	loader.RegisterFunction(range_fn);

	// ST_MaskBandValuesEq: inclusion — keep pixels equal to `value` only; nodata elsewhere.
	ScalarFunction eq_fn("ST_MaskBandValuesEq",
	                     {LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::DOUBLE, LogicalType::DOUBLE},
	                     LogicalType::BLOB, MaskBandValuesEqFunction, nullptr, nullptr, nullptr, nullptr,
	                     LogicalType(LogicalTypeId::INVALID), FunctionStability::CONSISTENT,
	                     FunctionNullHandling::SPECIAL_HANDLING);
	loader.RegisterFunction(eq_fn);
}

} // namespace duckdb
