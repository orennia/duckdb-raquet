#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <limits>

namespace duckdb {
namespace raquet {

// Band info structure (v0.3.0 format)
struct BandInfo {
    std::string name;
    std::string type;
    double nodata;
    bool has_nodata;

    BandInfo() : nodata(0), has_nodata(false) {}
    BandInfo(const std::string &n, const std::string &t)
        : name(n), type(t), nodata(0), has_nodata(false) {}
    BandInfo(const std::string &n, const std::string &t, double nd)
        : name(n), type(t), nodata(nd), has_nodata(true) {}
};

// Parsed raquet metadata (v0.3.0 format)
struct RaquetMetadata {
    std::string file_format;  // new in v0.3.0: should be "raquet"
    std::string compression;
    int block_width;
    int block_height;
    int min_zoom;       // was minresolution
    int max_zoom;       // was maxresolution
    int pixel_zoom;     // new in v0.3.0
    int num_blocks;
    std::string scheme; // "quadbin"
    std::string crs;    // "EPSG:3857"
    std::vector<BandInfo> band_info;  // Full band info including nodata
    std::vector<std::pair<std::string, std::string>> bands;  // name -> type (for backward compat)

    // Get band type by index (0-based) or name
    std::string get_band_type(int band_index) const {
        if (band_index < 0 || band_index >= static_cast<int>(bands.size())) {
            throw std::invalid_argument("Band index out of range");
        }
        return bands[band_index].second;
    }

    std::string get_band_type(const std::string &band_name) const {
        for (const auto &band : bands) {
            if (band.first == band_name) {
                return band.second;
            }
        }
        throw std::invalid_argument("Band not found: " + band_name);
    }

    // Get band info by index
    const BandInfo& get_band_info(int band_index) const {
        if (band_index < 0 || band_index >= static_cast<int>(band_info.size())) {
            throw std::invalid_argument("Band index out of range");
        }
        return band_info[band_index];
    }

    // Check if a value is nodata for a band
    bool is_nodata(int band_index, double value) const {
        if (band_index < 0 || band_index >= static_cast<int>(band_info.size())) {
            return false;
        }
        const auto &info = band_info[band_index];
        if (!info.has_nodata) return false;
        // Use exact comparison for integer types, approximate for float
        return value == info.nodata ||
               (std::isnan(value) && std::isnan(info.nodata));
    }
};

// Simple JSON value extraction (finds "key": value or "key": "value")
inline std::string extract_json_string(const std::string &json, const std::string &key) {
    std::string search = "\"" + key + "\":";
    size_t pos = json.find(search);
    if (pos == std::string::npos) {
        return "";
    }
    pos += search.length();

    // Skip whitespace
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t')) {
        pos++;
    }

    if (pos >= json.size()) return "";

    if (json[pos] == '"') {
        // String value
        pos++;
        size_t end = json.find('"', pos);
        if (end == std::string::npos) return "";
        return json.substr(pos, end - pos);
    } else {
        // Numeric or other value
        size_t end = pos;
        while (end < json.size() && json[end] != ',' && json[end] != '}' && json[end] != ']') {
            end++;
        }
        std::string val = json.substr(pos, end - pos);
        // Trim whitespace
        while (!val.empty() && (val.back() == ' ' || val.back() == '\t')) {
            val.pop_back();
        }
        return val;
    }
}

inline int extract_json_int(const std::string &json, const std::string &key, int default_val = 0) {
    std::string val = extract_json_string(json, key);
    if (val.empty()) return default_val;
    try {
        return std::stoi(val);
    } catch (...) {
        return default_val;
    }
}

inline double extract_json_double(const std::string &json, const std::string &key, double default_val = 0.0) {
    std::string val = extract_json_string(json, key);
    if (val.empty() || val == "null") return default_val;
    try {
        return std::stod(val);
    } catch (...) {
        return default_val;
    }
}

inline bool extract_json_has_value(const std::string &json, const std::string &key) {
    std::string val = extract_json_string(json, key);
    return !val.empty() && val != "null";
}

// Parse nodata value handling Zarr v3 string conventions:
// - "NaN" -> NaN
// - "Infinity" -> +Infinity
// - "-Infinity" -> -Infinity
// - numeric values as-is
inline double parse_nodata_value(const std::string &val) {
    if (val.empty() || val == "null") {
        return 0.0;
    }
    // Handle Zarr v3 string conventions (case-sensitive per spec)
    if (val == "NaN") {
        return std::nan("");
    }
    if (val == "Infinity") {
        return std::numeric_limits<double>::infinity();
    }
    if (val == "-Infinity") {
        return -std::numeric_limits<double>::infinity();
    }
    // Try to parse as numeric
    try {
        return std::stod(val);
    } catch (...) {
        return 0.0;
    }
}

// Extract a nested JSON object as a string
inline std::string extract_json_object(const std::string &json, const std::string &key) {
    std::string search = "\"" + key + "\":";
    size_t pos = json.find(search);
    if (pos == std::string::npos) return "";
    pos += search.length();

    // Skip whitespace
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t' || json[pos] == '\n')) {
        pos++;
    }

    if (pos >= json.size() || json[pos] != '{') return "";

    // Find matching closing brace
    int depth = 1;
    size_t start = pos;
    pos++;
    while (pos < json.size() && depth > 0) {
        if (json[pos] == '{') depth++;
        else if (json[pos] == '}') depth--;
        pos++;
    }

    return json.substr(start, pos - start);
}

// Parse bands array from metadata JSON (returns both legacy format and full BandInfo)
inline void parse_bands_full(const std::string &json,
                             std::vector<std::pair<std::string, std::string>> &bands,
                             std::vector<BandInfo> &band_info) {
    bands.clear();
    band_info.clear();

    // Find "bands": [...]
    size_t bands_pos = json.find("\"bands\":");
    if (bands_pos == std::string::npos) return;

    size_t arr_start = json.find('[', bands_pos);
    if (arr_start == std::string::npos) return;

    size_t arr_end = json.find(']', arr_start);
    if (arr_end == std::string::npos) return;

    std::string bands_str = json.substr(arr_start + 1, arr_end - arr_start - 1);

    // Parse each band object
    size_t pos = 0;
    while (pos < bands_str.size()) {
        size_t obj_start = bands_str.find('{', pos);
        if (obj_start == std::string::npos) break;

        size_t obj_end = bands_str.find('}', obj_start);
        if (obj_end == std::string::npos) break;

        std::string band_obj = bands_str.substr(obj_start, obj_end - obj_start + 1);

        std::string name = extract_json_string(band_obj, "name");
        std::string type = extract_json_string(band_obj, "type");

        if (!name.empty() && !type.empty()) {
            bands.push_back({name, type});

            BandInfo info(name, type);
            if (extract_json_has_value(band_obj, "nodata")) {
                // Use parse_nodata_value to handle Zarr v3 string conventions
                std::string nodata_str = extract_json_string(band_obj, "nodata");
                info.nodata = parse_nodata_value(nodata_str);
                info.has_nodata = true;
            }
            band_info.push_back(info);
        }

        pos = obj_end + 1;
    }
}

// Legacy function for backward compatibility
inline std::vector<std::pair<std::string, std::string>> parse_bands(const std::string &json) {
    std::vector<std::pair<std::string, std::string>> bands;
    std::vector<BandInfo> band_info;
    parse_bands_full(json, bands, band_info);
    return bands;
}

// Parse metadata JSON string (v0.3.0 format)
inline RaquetMetadata parse_metadata(const std::string &json) {
    RaquetMetadata meta;

    // New in v0.3.0: file_format identifier
    meta.file_format = extract_json_string(json, "file_format");

    meta.compression = extract_json_string(json, "compression");
    if (meta.compression.empty()) meta.compression = "none";

    meta.crs = extract_json_string(json, "crs");

    // Parse tiling object (v0.3.0)
    std::string tiling = extract_json_object(json, "tiling");
    if (!tiling.empty()) {
        meta.min_zoom = extract_json_int(tiling, "min_zoom", 0);
        meta.max_zoom = extract_json_int(tiling, "max_zoom", 26);
        meta.pixel_zoom = extract_json_int(tiling, "pixel_zoom", 0);
        meta.num_blocks = extract_json_int(tiling, "num_blocks", 0);
        meta.block_width = extract_json_int(tiling, "block_width", 256);
        meta.block_height = extract_json_int(tiling, "block_height", 256);
        meta.scheme = extract_json_string(tiling, "scheme");
    } else {
        // Defaults if no tiling object
        meta.min_zoom = 0;
        meta.max_zoom = 26;
        meta.pixel_zoom = 0;
        meta.num_blocks = 0;
        meta.block_width = 256;
        meta.block_height = 256;
        meta.scheme = "quadbin";
    }

    parse_bands_full(json, meta.bands, meta.band_info);

    return meta;
}

} // namespace raquet
} // namespace duckdb
