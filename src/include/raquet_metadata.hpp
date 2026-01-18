#pragma once

#include <string>
#include <vector>
#include <stdexcept>

namespace duckdb {
namespace raquet {

// Parsed raquet metadata
struct RaquetMetadata {
    std::string compression;
    int width;
    int height;
    int block_width;
    int block_height;
    int minresolution;
    int maxresolution;
    std::vector<std::pair<std::string, std::string>> bands;  // name -> type

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

// Parse bands array from metadata JSON
inline std::vector<std::pair<std::string, std::string>> parse_bands(const std::string &json) {
    std::vector<std::pair<std::string, std::string>> bands;

    // Find "bands": [...]
    size_t bands_pos = json.find("\"bands\":");
    if (bands_pos == std::string::npos) return bands;

    size_t arr_start = json.find('[', bands_pos);
    if (arr_start == std::string::npos) return bands;

    size_t arr_end = json.find(']', arr_start);
    if (arr_end == std::string::npos) return bands;

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
        }

        pos = obj_end + 1;
    }

    return bands;
}

// Parse metadata JSON string
inline RaquetMetadata parse_metadata(const std::string &json) {
    RaquetMetadata meta;

    meta.compression = extract_json_string(json, "compression");
    if (meta.compression.empty()) meta.compression = "none";

    meta.width = extract_json_int(json, "width", 256);
    meta.height = extract_json_int(json, "height", 256);
    meta.block_width = extract_json_int(json, "block_width", 256);
    meta.block_height = extract_json_int(json, "block_height", 256);
    meta.minresolution = extract_json_int(json, "minresolution", 0);
    meta.maxresolution = extract_json_int(json, "maxresolution", 26);

    meta.bands = parse_bands(json);

    return meta;
}

} // namespace raquet
} // namespace duckdb
