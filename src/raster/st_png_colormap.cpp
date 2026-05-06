#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "band_decoder.hpp"
#include "band_encoder.hpp"
#include "raquet_metadata.hpp"
#include "yyjson.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <unordered_map>
#include <string>
#include <vector>

using namespace duckdb_yyjson;

namespace duckdb {

namespace {

struct MapRgb {
	uint8_t r {0};
	uint8_t g {0};
	uint8_t b {0};
};

struct ColorMapEntry {
	bool is_exact {true};
	double v0 {0};
	double v1 {0};
	uint8_t r {0};
	uint8_t g {0};
	uint8_t b {0};
};

static void RequireSequentialGzipOrNone(const raquet::RaquetMetadata &meta, const char *fn) {
	if (meta.band_layout != "sequential") {
		throw InvalidInputException(
		    "%s: only band_layout='sequential' per-band tiles are supported (got '%s')", fn, meta.band_layout.c_str());
	}
	if (meta.compression != "gzip" && meta.compression != "none" && !meta.compression.empty()) {
		throw InvalidInputException(
		    "%s: only compression 'gzip' or 'none' is supported for sequential bands (got '%s')", fn,
		    meta.compression.c_str());
	}
}

static uint8_t DoubleToUint8Pixel(double p) {
	if (!std::isfinite(p)) {
		return 0;
	}
	long double x = std::llround(static_cast<long double>(p));
	if (x < 0) {
		x = 0;
	}
	if (x > 255) {
		x = 255;
	}
	return static_cast<uint8_t>(x);
}

static std::vector<double> DecodeBandBlob(const string_t &band, const raquet::RaquetMetadata &meta,
                                            const std::string &dtype_str) {
	const bool compressed = (meta.compression == "gzip");
	return raquet::decode_band(reinterpret_cast<const uint8_t *>(band.GetData()), band.GetSize(), dtype_str,
	                           meta.block_width, meta.block_height, compressed);
}

static string_t EncodeUint8TileBlob(Vector &blob_vec, const std::vector<uint8_t> &raw,
                                    const raquet::RaquetMetadata &meta) {
	std::vector<uint8_t> out;
	if (meta.compression == "gzip") {
		out = raquet::compress_gzip(raw.data(), raw.size());
	} else {
		out = raw;
	}
	return StringVector::AddStringOrBlob(blob_vec, reinterpret_cast<const char *>(out.data()), out.size());
}

static uint8_t JsonRgbComponent(yyjson_val *o, const char *key) {
	yyjson_val *v = yyjson_obj_get(o, key);
	if (!v || !yyjson_is_num(v)) {
		throw InvalidInputException("ST_ColorMapRaquet: entry or unmapped missing numeric '%s'", key);
	}
	const double d = yyjson_get_num(v);
	if (!std::isfinite(d) || d < 0 || d > 255) {
		throw InvalidInputException("ST_ColorMapRaquet: color component '%s' must be in [0,255]", key);
	}
	return static_cast<uint8_t>(std::llround(d));
}

static MapRgb ParseUnmapped(yyjson_val *root) {
	MapRgb out {0, 0, 0};
	yyjson_val *um = yyjson_obj_get(root, "unmapped");
	if (!um) {
		return out;
	}
	if (yyjson_is_str(um)) {
		const char *s = yyjson_get_str(um);
		if (std::strcmp(s, "nodata") == 0) {
			return out;
		}
		throw InvalidInputException(
		    "ST_ColorMapRaquet: unmapped string must be \"nodata\" or omit unmapped for default RGB 0,0,0");
	}
	if (!yyjson_is_obj(um)) {
		throw InvalidInputException("ST_ColorMapRaquet: unmapped must be \"nodata\" or an object {r,g,b}");
	}
	out.r = JsonRgbComponent(um, "r");
	out.g = JsonRgbComponent(um, "g");
	out.b = JsonRgbComponent(um, "b");
	return out;
}

static std::vector<ColorMapEntry> ParseColorMapEntriesRoot(yyjson_val *root) {
	if (!yyjson_is_obj(root)) {
		throw InvalidInputException("ST_ColorMapRaquet: colormap JSON must be an object");
	}
	yyjson_val *entries_val = yyjson_obj_get(root, "entries");
	if (!entries_val || !yyjson_is_arr(entries_val)) {
		throw InvalidInputException("ST_ColorMapRaquet: colormap must contain a non-empty \"entries\" array");
	}
	const size_t n = yyjson_arr_size(entries_val);
	if (n == 0) {
		throw InvalidInputException("ST_ColorMapRaquet: \"entries\" must be non-empty");
	}
	std::vector<ColorMapEntry> entries;
	entries.reserve(n);
	for (size_t i = 0; i < n; i++) {
		yyjson_val *entry = yyjson_arr_get(entries_val, i);
		if (!yyjson_is_obj(entry)) {
			throw InvalidInputException("ST_ColorMapRaquet: each entry must be an object");
		}
		ColorMapEntry e;
		e.r = JsonRgbComponent(entry, "r");
		e.g = JsonRgbComponent(entry, "g");
		e.b = JsonRgbComponent(entry, "b");
		yyjson_val *jv = yyjson_obj_get(entry, "value");
		yyjson_val *jlo = yyjson_obj_get(entry, "lo");
		yyjson_val *jhi = yyjson_obj_get(entry, "hi");
		const bool has_val = jv && yyjson_is_num(jv);
		const bool has_lo = jlo && yyjson_is_num(jlo);
		const bool has_hi = jhi && yyjson_is_num(jhi);
		if (has_val && (has_lo || has_hi)) {
			throw InvalidInputException("ST_ColorMapRaquet: entry cannot mix \"value\" with \"lo\"/\"hi\"");
		}
		if (has_val) {
			e.is_exact = true;
			e.v0 = yyjson_get_num(jv);
		} else if (has_lo && has_hi) {
			e.is_exact = false;
			e.v0 = yyjson_get_num(jlo);
			e.v1 = yyjson_get_num(jhi);
			if (!(e.v0 <= e.v1)) {
				throw InvalidInputException("ST_ColorMapRaquet: range entry requires lo <= hi");
			}
		} else {
			throw InvalidInputException("ST_ColorMapRaquet: each entry needs \"value\" or both \"lo\" and \"hi\"");
		}
		entries.push_back(e);
	}
	return entries;
}

static MapRgb MapPixel(double p, double nodata, const std::vector<ColorMapEntry> &entries, const MapRgb &unmapped) {
	if (!std::isfinite(p) || p == nodata) {
		return unmapped;
	}
	for (const auto &e : entries) {
		if (e.is_exact) {
			if (p == e.v0) {
				return MapRgb {e.r, e.g, e.b};
			}
		} else {
			if (p >= e.v0 && p <= e.v1) {
				return MapRgb {e.r, e.g, e.b};
			}
		}
	}
	return unmapped;
}

static int FindBandIndexByName(const raquet::RaquetMetadata &meta, const std::string &want_name) {
	for (size_t j = 0; j < meta.band_info.size(); j++) {
		if (meta.band_info[j].name == want_name) {
			return static_cast<int>(j);
		}
	}
	for (size_t j = 0; j < meta.bands.size(); j++) {
		if (meta.bands[j].first == want_name) {
			return static_cast<int>(j);
		}
	}
	return -1;
}

static MapRgb MapPixelPalette(double p,
                              const std::unordered_map<int64_t, std::array<uint8_t, 4>> &palette) {
	if (!std::isfinite(p)) {
		return MapRgb {0, 0, 0};
	}
	const int64_t k = static_cast<int64_t>(std::llround(p));
	auto it = palette.find(k);
	if (it == palette.end()) {
		return MapRgb {0, 0, 0};
	}
	const auto &rgba = it->second;
	return MapRgb {rgba[0], rgba[1], rgba[2]};
}

static void STAsPngRaquetFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	args.data[0].Flatten(args.size());
	args.data[1].Flatten(args.size());

	auto list_entries = ListVector::GetData(args.data[0]);
	auto &list_child = ListVector::GetEntry(args.data[0]);
	auto band_strings = FlatVector::GetData<string_t>(list_child);
	auto metadata_data = FlatVector::GetData<string_t>(args.data[1]);

	auto &list_validity = FlatVector::Validity(args.data[0]);
	auto &meta_validity = FlatVector::Validity(args.data[1]);

	auto result_data = FlatVector::GetData<string_t>(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!list_validity.RowIsValid(i) || !meta_validity.RowIsValid(i)) {
			FlatVector::Validity(result).SetInvalid(i);
			continue;
		}
		const auto &le = list_entries[i];
		if (le.length != 1 && le.length != 3 && le.length != 4) {
			throw InvalidInputException(
			    "ST_AsPNGRaquet: bands list must have length 1 (grayscale), 3 (RGB), or 4 (RGBA) (got %llu)",
			    static_cast<unsigned long long>(le.length));
		}
		raquet::RaquetMetadata meta;
		try {
			meta = raquet::parse_metadata(metadata_data[i].GetString());
		} catch (const std::exception &e) {
			throw InvalidInputException("ST_AsPNGRaquet: bad metadata: %s", e.what());
		}
		RequireSequentialGzipOrNone(meta, "ST_AsPNGRaquet");

		const int w = meta.block_width;
		const int h = meta.block_height;
		const size_t npix = static_cast<size_t>(w) * static_cast<size_t>(h);

		std::vector<std::vector<uint8_t>> plane_bytes;
		plane_bytes.reserve(le.length);
		bool row_ok = true;
		try {
			for (idx_t b = 0; b < le.length; b++) {
				const auto blob = band_strings[le.offset + b];
				if (blob.GetSize() == 0) {
					row_ok = false;
					break;
				}
				// RGB/RGBA lists are always uint8 tile planes (e.g. from colormap helpers); do not use
				// meta.bands[b] dtype, which may describe a single source band (e.g. float32) for all slots.
				const std::string dtype =
				    (le.length == 3 || le.length == 4)
				        ? std::string("uint8")
				        : (meta.bands.empty()
				               ? "uint8"
				               : ((static_cast<size_t>(b) < meta.bands.size()) ? meta.bands[b].second
				                                                              : meta.bands[0].second));
				auto pix = DecodeBandBlob(blob, meta, dtype);
				if (pix.size() != npix) {
					throw InvalidInputException("ST_AsPNGRaquet: decoded band size mismatch for band %llu",
					                            static_cast<unsigned long long>(b));
				}
				std::vector<uint8_t> raw(npix);
				for (size_t p = 0; p < npix; p++) {
					raw[p] = DoubleToUint8Pixel(pix[p]);
				}
				plane_bytes.push_back(std::move(raw));
			}
			if (!row_ok) {
				FlatVector::Validity(result).SetInvalid(i);
				continue;
			}
			std::vector<uint8_t> png_bytes;
			if (plane_bytes.size() == 1) {
				png_bytes = raquet::encode_png(plane_bytes[0].data(), w, h, 1);
			} else if (plane_bytes.size() == 3) {
				auto inter = raquet::interleave_bands(plane_bytes, w, h, 1);
				png_bytes = raquet::encode_png(inter.data(), w, h, 3);
			} else {
				auto inter = raquet::interleave_bands(plane_bytes, w, h, 1);
				png_bytes = raquet::encode_png(inter.data(), w, h, 4);
			}
			result_data[i] =
			    StringVector::AddStringOrBlob(result, reinterpret_cast<const char *>(png_bytes.data()), png_bytes.size());
		} catch (const InvalidInputException &) {
			throw;
		} catch (...) {
			FlatVector::Validity(result).SetInvalid(i);
		}
	}
}

static void STColorMapRaquetFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	args.data[0].Flatten(args.size());
	args.data[1].Flatten(args.size());
	args.data[2].Flatten(args.size());
	args.data[3].Flatten(args.size());

	auto band_data = FlatVector::GetData<string_t>(args.data[0]);
	auto metadata_data = FlatVector::GetData<string_t>(args.data[1]);
	auto nodata_data = FlatVector::GetData<double>(args.data[2]);
	auto map_data = FlatVector::GetData<string_t>(args.data[3]);

	auto &band_validity = FlatVector::Validity(args.data[0]);
	auto &meta_validity = FlatVector::Validity(args.data[1]);
	auto &nodata_validity = FlatVector::Validity(args.data[2]);
	auto &map_validity = FlatVector::Validity(args.data[3]);

	auto &struct_entries = StructVector::GetEntries(result);
	auto r_out = FlatVector::GetData<string_t>(*struct_entries[0]);
	auto g_out = FlatVector::GetData<string_t>(*struct_entries[1]);
	auto b_out = FlatVector::GetData<string_t>(*struct_entries[2]);
	auto &result_validity = FlatVector::Validity(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!band_validity.RowIsValid(i) || !meta_validity.RowIsValid(i) || !nodata_validity.RowIsValid(i) ||
		    !map_validity.RowIsValid(i)) {
			result_validity.SetInvalid(i);
			continue;
		}
		const auto band = band_data[i];
		if (band.GetSize() == 0) {
			result_validity.SetInvalid(i);
			continue;
		}
		raquet::RaquetMetadata meta;
		try {
			meta = raquet::parse_metadata(metadata_data[i].GetString());
		} catch (const std::exception &e) {
			throw InvalidInputException("ST_ColorMapRaquet: bad metadata: %s", e.what());
		}
		RequireSequentialGzipOrNone(meta, "ST_ColorMapRaquet");

		const double nodata = nodata_data[i];
		const std::string map_str = map_data[i].GetString();
		yyjson_doc *doc = yyjson_read(map_str.c_str(), map_str.length(), 0);
		if (!doc) {
			throw InvalidInputException("ST_ColorMapRaquet: invalid colormap JSON");
		}
		yyjson_val *root = yyjson_doc_get_root(doc);
		MapRgb unmapped {0, 0, 0};
		std::vector<ColorMapEntry> entries;
		try {
			if (!yyjson_is_obj(root)) {
				throw InvalidInputException("ST_ColorMapRaquet: colormap JSON must be an object");
			}
			unmapped = ParseUnmapped(root);
			entries = ParseColorMapEntriesRoot(root);
		} catch (...) {
			yyjson_doc_free(doc);
			throw;
		}
		yyjson_doc_free(doc);

		const std::string dtype0 = meta.bands.empty() ? std::string("uint8") : meta.bands[0].second;
		std::vector<double> pixels;
		try {
			pixels = DecodeBandBlob(band, meta, dtype0);
		} catch (...) {
			result_validity.SetInvalid(i);
			continue;
		}
		const size_t npix = pixels.size();
		std::vector<uint8_t> rr(npix), gg(npix), bb(npix);
		for (size_t p = 0; p < npix; p++) {
			const MapRgb c = MapPixel(pixels[p], nodata, entries, unmapped);
			rr[p] = c.r;
			gg[p] = c.g;
			bb[p] = c.b;
		}
		try {
			r_out[i] = EncodeUint8TileBlob(*struct_entries[0], rr, meta);
			g_out[i] = EncodeUint8TileBlob(*struct_entries[1], gg, meta);
			b_out[i] = EncodeUint8TileBlob(*struct_entries[2], bb, meta);
		} catch (...) {
			result_validity.SetInvalid(i);
		}
	}

	result.SetVectorType(VectorType::FLAT_VECTOR);
}

static void STColorMapRaquetMetadataOnlyFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	args.data[0].Flatten(args.size());
	args.data[1].Flatten(args.size());

	auto band_data = FlatVector::GetData<string_t>(args.data[0]);
	auto metadata_data = FlatVector::GetData<string_t>(args.data[1]);

	auto &band_validity = FlatVector::Validity(args.data[0]);
	auto &meta_validity = FlatVector::Validity(args.data[1]);

	auto &struct_entries = StructVector::GetEntries(result);
	auto r_out = FlatVector::GetData<string_t>(*struct_entries[0]);
	auto g_out = FlatVector::GetData<string_t>(*struct_entries[1]);
	auto b_out = FlatVector::GetData<string_t>(*struct_entries[2]);
	auto &result_validity = FlatVector::Validity(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!band_validity.RowIsValid(i) || !meta_validity.RowIsValid(i)) {
			result_validity.SetInvalid(i);
			continue;
		}
		const auto band = band_data[i];
		if (band.GetSize() == 0) {
			result_validity.SetInvalid(i);
			continue;
		}
		raquet::RaquetMetadata meta;
		try {
			meta = raquet::parse_metadata(metadata_data[i].GetString());
		} catch (const std::exception &e) {
			throw InvalidInputException("ST_ColorMapRaquet: bad metadata: %s", e.what());
		}
		RequireSequentialGzipOrNone(meta, "ST_ColorMapRaquet");

		if (meta.band_info.empty() || !meta.band_info[0].has_colortable ||
		    meta.band_info[0].colortable.empty()) {
			throw InvalidInputException(
			    "ST_ColorMapRaquet(band, metadata): first band metadata must contain a non-empty \"colortable\" object — "
			    "otherwise use ST_ColorMapRaquet(band, metadata, nodata, map_json) or ST_ColorMapRaquet(band, band_name, metadata)");
		}

		const auto &palette = meta.band_info[0].colortable;
		const std::string dtype0 =
		    meta.band_info[0].type.empty()
		        ? (meta.bands.empty() ? std::string("uint8") : meta.bands[0].second)
		        : meta.band_info[0].type;
		std::vector<double> pixels;
		try {
			pixels = DecodeBandBlob(band, meta, dtype0);
		} catch (...) {
			result_validity.SetInvalid(i);
			continue;
		}

		const size_t npix = pixels.size();
		std::vector<uint8_t> rr(npix), gg(npix), bb(npix);
		for (size_t p = 0; p < npix; p++) {
			const MapRgb c = MapPixelPalette(pixels[p], palette);
			rr[p] = c.r;
			gg[p] = c.g;
			bb[p] = c.b;
		}
		try {
			r_out[i] = EncodeUint8TileBlob(*struct_entries[0], rr, meta);
			g_out[i] = EncodeUint8TileBlob(*struct_entries[1], gg, meta);
			b_out[i] = EncodeUint8TileBlob(*struct_entries[2], bb, meta);
		} catch (...) {
			result_validity.SetInvalid(i);
		}
	}

	result.SetVectorType(VectorType::FLAT_VECTOR);
}

static void STColorMapRaquetNamedBandFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	args.data[0].Flatten(args.size());
	args.data[1].Flatten(args.size());
	args.data[2].Flatten(args.size());

	auto band_data = FlatVector::GetData<string_t>(args.data[0]);
	auto band_name_data = FlatVector::GetData<string_t>(args.data[1]);
	auto metadata_data = FlatVector::GetData<string_t>(args.data[2]);

	auto &band_validity = FlatVector::Validity(args.data[0]);
	auto &name_validity = FlatVector::Validity(args.data[1]);
	auto &meta_validity = FlatVector::Validity(args.data[2]);

	auto &struct_entries = StructVector::GetEntries(result);
	auto r_out = FlatVector::GetData<string_t>(*struct_entries[0]);
	auto g_out = FlatVector::GetData<string_t>(*struct_entries[1]);
	auto b_out = FlatVector::GetData<string_t>(*struct_entries[2]);
	auto &result_validity = FlatVector::Validity(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!band_validity.RowIsValid(i) || !name_validity.RowIsValid(i) || !meta_validity.RowIsValid(i)) {
			result_validity.SetInvalid(i);
			continue;
		}
		const auto band = band_data[i];
		if (band.GetSize() == 0) {
			result_validity.SetInvalid(i);
			continue;
		}
		const std::string want_name = band_name_data[i].GetString();
		if (want_name.empty()) {
			throw InvalidInputException("ST_ColorMapRaquet(band, band_name, metadata): band_name must be non-empty");
		}

		raquet::RaquetMetadata meta;
		try {
			meta = raquet::parse_metadata(metadata_data[i].GetString());
		} catch (const std::exception &e) {
			throw InvalidInputException("ST_ColorMapRaquet: bad metadata: %s", e.what());
		}
		RequireSequentialGzipOrNone(meta, "ST_ColorMapRaquet");

		const int bi = FindBandIndexByName(meta, want_name);
		if (bi < 0) {
			throw InvalidInputException(
			    "ST_ColorMapRaquet(band, band_name, metadata): no band named '%s' in metadata "
			    "(use the same name as in bands[].name)",
			    want_name.c_str());
		}
		if (bi >= static_cast<int>(meta.band_info.size())) {
			throw InvalidInputException(
			    "ST_ColorMapRaquet(band, band_name, metadata): band '%s' has no full band object in metadata",
			    want_name.c_str());
		}
		const auto &binfo = meta.band_info[static_cast<size_t>(bi)];
		if (!binfo.has_colortable || binfo.colortable.empty()) {
			throw InvalidInputException(
			    "ST_ColorMapRaquet(band, band_name, metadata): band '%s' has no non-empty \"colortable\" in metadata",
			    want_name.c_str());
		}

		const auto &palette = binfo.colortable;
		const std::string dtype_str =
		    binfo.type.empty() ? (static_cast<size_t>(bi) < meta.bands.size() ? meta.bands[static_cast<size_t>(bi)].second
		                                                                          : std::string("uint8"))
		                       : binfo.type;

		std::vector<double> pixels;
		try {
			pixels = DecodeBandBlob(band, meta, dtype_str);
		} catch (...) {
			result_validity.SetInvalid(i);
			continue;
		}

		const size_t npix = pixels.size();
		std::vector<uint8_t> rr(npix), gg(npix), bb(npix);
		for (size_t p = 0; p < npix; p++) {
			const MapRgb c = MapPixelPalette(pixels[p], palette);
			rr[p] = c.r;
			gg[p] = c.g;
			bb[p] = c.b;
		}
		try {
			r_out[i] = EncodeUint8TileBlob(*struct_entries[0], rr, meta);
			g_out[i] = EncodeUint8TileBlob(*struct_entries[1], gg, meta);
			b_out[i] = EncodeUint8TileBlob(*struct_entries[2], bb, meta);
		} catch (...) {
			result_validity.SetInvalid(i);
		}
	}

	result.SetVectorType(VectorType::FLAT_VECTOR);
}

// --- ST_ColorMapContinuousRaquet: piecewise-linear ramp in RGB (stops in data value space) ---

struct ContinuousSpec {
	std::vector<std::pair<double, MapRgb>> stops;
	bool has_vmin {false};
	bool has_vmax {false};
	double vmin {0};
	double vmax {1};
	bool has_effective_nodata {false};
	double effective_nodata {0};
	bool outside_clamp {true};
	MapRgb unmapped {0, 0, 0};
};

static uint8_t ParseHexByte(const char *s, size_t &i, size_t n) {
	if (i + 1 >= n) {
		throw InvalidInputException("ST_ColorMapContinuousRaquet: invalid hex color");
	}
	auto hexval = [](char c) -> int {
		if (c >= '0' && c <= '9') {
			return c - '0';
		}
		if (c >= 'a' && c <= 'f') {
			return 10 + c - 'a';
		}
		if (c >= 'A' && c <= 'F') {
			return 10 + c - 'A';
		}
		throw InvalidInputException("ST_ColorMapContinuousRaquet: invalid hex digit");
	};
	const int hi = hexval(s[i++]);
	const int lo = hexval(s[i++]);
	return static_cast<uint8_t>(hi * 16 + lo);
}

static MapRgb ParseHexColorStr(const char *str, size_t len) {
	if (!str || len < 2 || str[0] != '#') {
		throw InvalidInputException("ST_ColorMapContinuousRaquet: color hex must start with '#'");
	}
	MapRgb out {0, 0, 0};
	auto hex1 = [](char c) -> int {
		if (c >= '0' && c <= '9') {
			return c - '0';
		}
		if (c >= 'a' && c <= 'f') {
			return 10 + c - 'a';
		}
		if (c >= 'A' && c <= 'F') {
			return 10 + c - 'A';
		}
		throw InvalidInputException("ST_ColorMapContinuousRaquet: invalid hex digit in color");
	};
	if (len == 4) {
		// #RGB
		out.r = static_cast<uint8_t>(hex1(str[1]) * 17);
		out.g = static_cast<uint8_t>(hex1(str[2]) * 17);
		out.b = static_cast<uint8_t>(hex1(str[3]) * 17);
		return out;
	}
	if (len >= 7) {
		// #RRGGBB or #RRGGBBAA (alpha ignored for RGB output)
		size_t i = 1;
		out.r = ParseHexByte(str, i, len);
		out.g = ParseHexByte(str, i, len);
		out.b = ParseHexByte(str, i, len);
		return out;
	}
	throw InvalidInputException("ST_ColorMapContinuousRaquet: unsupported hex color length");
}

static MapRgb ParseStopColor(yyjson_val *obj, const char *ctx) {
	yyjson_val *hexv = yyjson_obj_get(obj, "hex");
	if (hexv && yyjson_is_str(hexv)) {
		const char *hs = yyjson_get_str(hexv);
		const size_t ln = std::strlen(hs);
		return ParseHexColorStr(hs, ln);
	}
	yyjson_val *jr = yyjson_obj_get(obj, "r");
	yyjson_val *jg = yyjson_obj_get(obj, "g");
	yyjson_val *jb = yyjson_obj_get(obj, "b");
	if (jr && jg && jb && yyjson_is_num(jr) && yyjson_is_num(jg) && yyjson_is_num(jb)) {
		const double dr = yyjson_get_num(jr);
		const double dg = yyjson_get_num(jg);
		const double db = yyjson_get_num(jb);
		if (!std::isfinite(dr) || !std::isfinite(dg) || !std::isfinite(db) || dr < 0 || dr > 255 || dg < 0 || dg > 255 ||
		    db < 0 || db > 255) {
			throw InvalidInputException("%s: stop r,g,b must be finite and in [0,255]", ctx);
		}
		return MapRgb {static_cast<uint8_t>(std::llround(dr)), static_cast<uint8_t>(std::llround(dg)),
		               static_cast<uint8_t>(std::llround(db))};
	}
	throw InvalidInputException(
	    "%s: each stop needs \"hex\" (e.g. \"#RRGGBB\") or numeric \"r\",\"g\",\"b\" in [0,255]", ctx);
}

static double JsonStopV(yyjson_val *obj) {
	yyjson_val *jv = yyjson_obj_get(obj, "v");
	if (!jv) {
		jv = yyjson_obj_get(obj, "value");
	}
	if (!jv) {
		jv = yyjson_obj_get(obj, "at");
	}
	if (!jv || !yyjson_is_num(jv)) {
		throw InvalidInputException(
		    "ST_ColorMapContinuousRaquet: each stop must have a numeric \"v\", \"value\", or \"at\"");
	}
	const double v = yyjson_get_num(jv);
	if (!std::isfinite(v)) {
		throw InvalidInputException("ST_ColorMapContinuousRaquet: stop position must be finite");
	}
	return v;
}

static ContinuousSpec ParseContinuousSpecRoot(yyjson_val *root, const raquet::BandInfo *band_for_nodata) {
	ContinuousSpec spec;
	if (!yyjson_is_obj(root)) {
		throw InvalidInputException("ST_ColorMapContinuousRaquet: spec must be a JSON object");
	}
	yyjson_val *stops_arr = yyjson_obj_get(root, "stops");
	if (!stops_arr || !yyjson_is_arr(stops_arr) || yyjson_arr_size(stops_arr) < 2) {
		throw InvalidInputException("ST_ColorMapContinuousRaquet: spec must contain \"stops\" array with at least 2 entries");
	}
	const size_t ns = yyjson_arr_size(stops_arr);
	spec.stops.reserve(ns);
	for (size_t i = 0; i < ns; i++) {
		yyjson_val *el = yyjson_arr_get(stops_arr, i);
		if (!yyjson_is_obj(el)) {
			throw InvalidInputException("ST_ColorMapContinuousRaquet: each stop must be an object");
		}
		const double v = JsonStopV(el);
		const MapRgb c = ParseStopColor(el, "ST_ColorMapContinuousRaquet");
		spec.stops.emplace_back(v, c);
	}
	std::sort(spec.stops.begin(), spec.stops.end(),
	           [](const std::pair<double, MapRgb> &a, const std::pair<double, MapRgb> &b) { return a.first < b.first; });

	yyjson_val *jvm = yyjson_obj_get(root, "vmin");
	if (jvm && yyjson_is_num(jvm)) {
		spec.has_vmin = true;
		spec.vmin = yyjson_get_num(jvm);
	}
	yyjson_val *jvx = yyjson_obj_get(root, "vmax");
	if (jvx && yyjson_is_num(jvx)) {
		spec.has_vmax = true;
		spec.vmax = yyjson_get_num(jvx);
	}

	yyjson_val *jnd = yyjson_obj_get(root, "nodata");
	if (jnd && yyjson_is_null(jnd)) {
		// explicit: do not treat any value as nodata for this ramp
		spec.has_effective_nodata = false;
	} else if (jnd && yyjson_is_num(jnd)) {
		spec.has_effective_nodata = true;
		spec.effective_nodata = yyjson_get_num(jnd);
	} else if (band_for_nodata && band_for_nodata->has_nodata) {
		spec.has_effective_nodata = true;
		spec.effective_nodata = band_for_nodata->nodata;
	}

	yyjson_val *outv = yyjson_obj_get(root, "outside");
	if (outv && yyjson_is_str(outv)) {
		const char *os = yyjson_get_str(outv);
		if (std::strcmp(os, "nodata") == 0) {
			spec.outside_clamp = false;
		} else if (std::strcmp(os, "clamp") == 0) {
			spec.outside_clamp = true;
		} else {
			throw InvalidInputException("ST_ColorMapContinuousRaquet: \"outside\" must be \"clamp\" or \"nodata\"");
		}
	}

	yyjson_val *um = yyjson_obj_get(root, "unmapped");
	if (um && yyjson_is_obj(um)) {
		yyjson_val *jr = yyjson_obj_get(um, "r");
		yyjson_val *jg = yyjson_obj_get(um, "g");
		yyjson_val *jb = yyjson_obj_get(um, "b");
		if (!jr || !jg || !jb || !yyjson_is_num(jr) || !yyjson_is_num(jg) || !yyjson_is_num(jb)) {
			throw InvalidInputException("ST_ColorMapContinuousRaquet: \"unmapped\" must include numeric r, g, b");
		}
		const double dr = yyjson_get_num(jr);
		const double dg = yyjson_get_num(jg);
		const double db = yyjson_get_num(jb);
		if (!std::isfinite(dr) || !std::isfinite(dg) || !std::isfinite(db) || dr < 0 || dr > 255 || dg < 0 || dg > 255 ||
		    db < 0 || db > 255) {
			throw InvalidInputException("ST_ColorMapContinuousRaquet: unmapped r,g,b must be finite and in [0,255]");
		}
		spec.unmapped.r = static_cast<uint8_t>(std::llround(dr));
		spec.unmapped.g = static_cast<uint8_t>(std::llround(dg));
		spec.unmapped.b = static_cast<uint8_t>(std::llround(db));
	}

	return spec;
}

static void TileMinMaxFiniteExcludingNodata(const std::vector<double> &pixels, bool has_nd, double nd, double &out_min,
                                            double &out_max, bool &any) {
	any = false;
	out_min = std::numeric_limits<double>::infinity();
	out_max = -std::numeric_limits<double>::infinity();
	for (double p : pixels) {
		if (!std::isfinite(p)) {
			continue;
		}
		if (has_nd && p == nd) {
			continue;
		}
		any = true;
		out_min = std::min(out_min, p);
		out_max = std::max(out_max, p);
	}
}

static MapRgb LerpRgb(const MapRgb &a, const MapRgb &b, double t) {
	auto ch = [](uint8_t x, uint8_t y, double tt) {
		return static_cast<uint8_t>(
		    std::llround(static_cast<double>(x) + tt * (static_cast<double>(y) - static_cast<double>(x))));
	};
	return MapRgb {ch(a.r, b.r, t), ch(a.g, b.g, t), ch(a.b, b.b, t)};
}

static MapRgb ColorAlongStops(double q, const std::vector<std::pair<double, MapRgb>> &st) {
	if (st.empty()) {
		return MapRgb {0, 0, 0};
	}
	if (q <= st[0].first) {
		return st[0].second;
	}
	if (q >= st.back().first) {
		return st.back().second;
	}
	for (size_t i = 0; i + 1 < st.size(); i++) {
		const double v0 = st[i].first;
		const double v1 = st[i + 1].first;
		if (q >= v0 && q <= v1) {
			if (v1 == v0) {
				return st[i].second;
			}
			const double t = (q - v0) / (v1 - v0);
			return LerpRgb(st[i].second, st[i + 1].second, t);
		}
	}
	return st.back().second;
}

static MapRgb MapPixelContinuous(double p, const ContinuousSpec &spec, double eff_vmin, double eff_vmax) {
	if (!std::isfinite(p)) {
		return spec.unmapped;
	}
	if (spec.has_effective_nodata && p == spec.effective_nodata) {
		return spec.unmapped;
	}
	if (!spec.outside_clamp && (p < eff_vmin || p > eff_vmax)) {
		return spec.unmapped;
	}
	double q = p;
	if (spec.outside_clamp) {
		q = std::min(std::max(p, eff_vmin), eff_vmax);
	}
	if (!(eff_vmax > eff_vmin)) {
		return spec.stops[0].second;
	}
	return ColorAlongStops(q, spec.stops);
}

static void ApplyContinuousColormapRow(const string_t &band, const raquet::RaquetMetadata &meta, int band_ix,
                                       const std::string &spec_json, Vector &r_vec, Vector &g_vec, Vector &b_vec, idx_t row,
                                       string_t *r_out, string_t *g_out, string_t *b_out) {
	RequireSequentialGzipOrNone(meta, "ST_ColorMapContinuousRaquet");

	const raquet::BandInfo *binfo =
	    (band_ix >= 0 && band_ix < static_cast<int>(meta.band_info.size())) ? &meta.band_info[static_cast<size_t>(band_ix)]
	                                                                         : nullptr;
	std::string dtype_str = "float32";
	if (binfo && !binfo->type.empty()) {
		dtype_str = binfo->type;
	} else if (band_ix >= 0 && static_cast<size_t>(band_ix) < meta.bands.size()) {
		dtype_str = meta.bands[static_cast<size_t>(band_ix)].second;
	} else if (!meta.bands.empty()) {
		dtype_str = meta.bands[0].second;
	}

	yyjson_doc *doc = yyjson_read(spec_json.c_str(), spec_json.length(), 0);
	if (!doc) {
		throw InvalidInputException("ST_ColorMapContinuousRaquet: invalid spec JSON");
	}
	ContinuousSpec spec;
	try {
		spec = ParseContinuousSpecRoot(yyjson_doc_get_root(doc), binfo);
	} catch (...) {
		yyjson_doc_free(doc);
		throw;
	}
	yyjson_doc_free(doc);

	std::vector<double> pixels = DecodeBandBlob(band, meta, dtype_str);
	double tile_min = 0, tile_max = 1;
	bool any_finite = false;
	TileMinMaxFiniteExcludingNodata(pixels, spec.has_effective_nodata, spec.effective_nodata, tile_min, tile_max,
	                                any_finite);
	if (!any_finite) {
		tile_min = 0;
		tile_max = 1;
	}
	double eff_vmin = spec.has_vmin ? spec.vmin : tile_min;
	double eff_vmax = spec.has_vmax ? spec.vmax : tile_max;
	if (eff_vmin > eff_vmax) {
		std::swap(eff_vmin, eff_vmax);
	}

	const size_t npix = pixels.size();
	std::vector<uint8_t> rr(npix), gg(npix), bb(npix);
	for (size_t p = 0; p < npix; p++) {
		const MapRgb c = MapPixelContinuous(pixels[p], spec, eff_vmin, eff_vmax);
		rr[p] = c.r;
		gg[p] = c.g;
		bb[p] = c.b;
	}
	r_out[row] = EncodeUint8TileBlob(r_vec, rr, meta);
	g_out[row] = EncodeUint8TileBlob(g_vec, gg, meta);
	b_out[row] = EncodeUint8TileBlob(b_vec, bb, meta);
}

static void STColorMapContinuousRaquetFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	args.data[0].Flatten(args.size());
	args.data[1].Flatten(args.size());
	args.data[2].Flatten(args.size());

	auto band_data = FlatVector::GetData<string_t>(args.data[0]);
	auto metadata_data = FlatVector::GetData<string_t>(args.data[1]);
	auto spec_data = FlatVector::GetData<string_t>(args.data[2]);

	auto &band_validity = FlatVector::Validity(args.data[0]);
	auto &meta_validity = FlatVector::Validity(args.data[1]);
	auto &spec_validity = FlatVector::Validity(args.data[2]);

	auto &struct_entries = StructVector::GetEntries(result);
	auto r_out = FlatVector::GetData<string_t>(*struct_entries[0]);
	auto g_out = FlatVector::GetData<string_t>(*struct_entries[1]);
	auto b_out = FlatVector::GetData<string_t>(*struct_entries[2]);
	auto &result_validity = FlatVector::Validity(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!band_validity.RowIsValid(i) || !meta_validity.RowIsValid(i) || !spec_validity.RowIsValid(i)) {
			result_validity.SetInvalid(i);
			continue;
		}
		const auto band = band_data[i];
		if (band.GetSize() == 0) {
			result_validity.SetInvalid(i);
			continue;
		}
		raquet::RaquetMetadata meta;
		try {
			meta = raquet::parse_metadata(metadata_data[i].GetString());
		} catch (const std::exception &e) {
			throw InvalidInputException("ST_ColorMapContinuousRaquet: bad metadata: %s", e.what());
		}
		try {
			ApplyContinuousColormapRow(band, meta, 0, spec_data[i].GetString(), *struct_entries[0], *struct_entries[1],
			                           *struct_entries[2], i, r_out, g_out, b_out);
		} catch (const InvalidInputException &) {
			throw;
		} catch (const std::exception &e) {
			throw InvalidInputException("ST_ColorMapContinuousRaquet: %s", e.what());
		}
	}
	result.SetVectorType(VectorType::FLAT_VECTOR);
}

static void STColorMapContinuousNamedRaquetFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	args.data[0].Flatten(args.size());
	args.data[1].Flatten(args.size());
	args.data[2].Flatten(args.size());
	args.data[3].Flatten(args.size());

	auto band_data = FlatVector::GetData<string_t>(args.data[0]);
	auto name_data = FlatVector::GetData<string_t>(args.data[1]);
	auto metadata_data = FlatVector::GetData<string_t>(args.data[2]);
	auto spec_data = FlatVector::GetData<string_t>(args.data[3]);

	auto &band_validity = FlatVector::Validity(args.data[0]);
	auto &name_validity = FlatVector::Validity(args.data[1]);
	auto &meta_validity = FlatVector::Validity(args.data[2]);
	auto &spec_validity = FlatVector::Validity(args.data[3]);

	auto &struct_entries = StructVector::GetEntries(result);
	auto r_out = FlatVector::GetData<string_t>(*struct_entries[0]);
	auto g_out = FlatVector::GetData<string_t>(*struct_entries[1]);
	auto b_out = FlatVector::GetData<string_t>(*struct_entries[2]);
	auto &result_validity = FlatVector::Validity(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!band_validity.RowIsValid(i) || !name_validity.RowIsValid(i) || !meta_validity.RowIsValid(i) ||
		    !spec_validity.RowIsValid(i)) {
			result_validity.SetInvalid(i);
			continue;
		}
		const auto band = band_data[i];
		if (band.GetSize() == 0) {
			result_validity.SetInvalid(i);
			continue;
		}
		const std::string want = name_data[i].GetString();
		if (want.empty()) {
			throw InvalidInputException("ST_ColorMapContinuousRaquet: band_name must be non-empty");
		}
		raquet::RaquetMetadata meta;
		try {
			meta = raquet::parse_metadata(metadata_data[i].GetString());
		} catch (const std::exception &e) {
			throw InvalidInputException("ST_ColorMapContinuousRaquet: bad metadata: %s", e.what());
		}
		const int bi = FindBandIndexByName(meta, want);
		if (bi < 0) {
			throw InvalidInputException("ST_ColorMapContinuousRaquet: unknown band name '%s'", want.c_str());
		}
		try {
			ApplyContinuousColormapRow(band, meta, bi, spec_data[i].GetString(), *struct_entries[0], *struct_entries[1],
			                           *struct_entries[2], i, r_out, g_out, b_out);
		} catch (const InvalidInputException &) {
			throw;
		} catch (const std::exception &e) {
			throw InvalidInputException("ST_ColorMapContinuousRaquet: %s", e.what());
		}
	}
	result.SetVectorType(VectorType::FLAT_VECTOR);
}

} // namespace

void RegisterPngColormapFunctions(ExtensionLoader &loader) {
	child_list_t<LogicalType> cmap_fields;
	cmap_fields.push_back(make_pair("band_r", LogicalType::BLOB));
	cmap_fields.push_back(make_pair("band_g", LogicalType::BLOB));
	cmap_fields.push_back(make_pair("band_b", LogicalType::BLOB));
	const auto cmap_struct = LogicalType::STRUCT(cmap_fields);

	ScalarFunction as_png_fn("ST_AsPNGRaquet",
	                         {LogicalType::LIST(LogicalType::BLOB), LogicalType::VARCHAR}, LogicalType::BLOB,
	                         STAsPngRaquetFunction, nullptr, nullptr, nullptr, nullptr, LogicalType(LogicalTypeId::INVALID),
	                         FunctionStability::CONSISTENT, FunctionNullHandling::DEFAULT_NULL_HANDLING);
	loader.RegisterFunction(as_png_fn);

	ScalarFunction cmap_meta_fn(
	    "ST_ColorMapRaquet", {LogicalType::BLOB, LogicalType::VARCHAR}, cmap_struct, STColorMapRaquetMetadataOnlyFunction,
	    nullptr, nullptr, nullptr, nullptr, LogicalType(LogicalTypeId::INVALID), FunctionStability::CONSISTENT,
	    FunctionNullHandling::DEFAULT_NULL_HANDLING);
	loader.RegisterFunction(cmap_meta_fn);

	ScalarFunction cmap_named_fn("ST_ColorMapRaquet",
	                             {LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::VARCHAR}, cmap_struct,
	                             STColorMapRaquetNamedBandFunction, nullptr, nullptr, nullptr, nullptr,
	                             LogicalType(LogicalTypeId::INVALID), FunctionStability::CONSISTENT,
	                             FunctionNullHandling::DEFAULT_NULL_HANDLING);
	loader.RegisterFunction(cmap_named_fn);

	ScalarFunction cmap_fn("ST_ColorMapRaquet",
	                       {LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::DOUBLE, LogicalType::VARCHAR},
	                       cmap_struct, STColorMapRaquetFunction, nullptr, nullptr, nullptr, nullptr,
	                       LogicalType(LogicalTypeId::INVALID), FunctionStability::CONSISTENT,
	                       FunctionNullHandling::DEFAULT_NULL_HANDLING);
	loader.RegisterFunction(cmap_fn);

	ScalarFunction cont_fn("ST_ColorMapContinuousRaquet",
	                       {LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::VARCHAR}, cmap_struct,
	                       STColorMapContinuousRaquetFunction, nullptr, nullptr, nullptr, nullptr,
	                       LogicalType(LogicalTypeId::INVALID), FunctionStability::CONSISTENT,
	                       FunctionNullHandling::DEFAULT_NULL_HANDLING);
	loader.RegisterFunction(cont_fn);

	ScalarFunction cont_named_fn(
	    "ST_ColorMapContinuousRaquet",
	    {LogicalType::BLOB, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}, cmap_struct,
	    STColorMapContinuousNamedRaquetFunction, nullptr, nullptr, nullptr, nullptr, LogicalType(LogicalTypeId::INVALID),
	    FunctionStability::CONSISTENT, FunctionNullHandling::DEFAULT_NULL_HANDLING);
	loader.RegisterFunction(cont_named_fn);
}

} // namespace duckdb
