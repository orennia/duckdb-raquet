/// Raquet metadata parsing.
/// Ported from the C++ raquet_metadata.hpp header.

/// Band data type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BandDataType {
    Uint8,
    Int8,
    Uint16,
    Int16,
    Uint32,
    Int32,
    Uint64,
    Int64,
    Float16,
    Float32,
    Float64,
}

impl BandDataType {
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "uint8" => Some(Self::Uint8),
            "int8" => Some(Self::Int8),
            "uint16" => Some(Self::Uint16),
            "int16" => Some(Self::Int16),
            "uint32" => Some(Self::Uint32),
            "int32" => Some(Self::Int32),
            "uint64" => Some(Self::Uint64),
            "int64" => Some(Self::Int64),
            "float16" => Some(Self::Float16),
            "float32" => Some(Self::Float32),
            "float64" => Some(Self::Float64),
            _ => None,
        }
    }

    pub fn byte_size(self) -> usize {
        match self {
            Self::Uint8 | Self::Int8 => 1,
            Self::Uint16 | Self::Int16 | Self::Float16 => 2,
            Self::Uint32 | Self::Int32 | Self::Float32 => 4,
            Self::Uint64 | Self::Int64 | Self::Float64 => 8,
        }
    }
}

/// Per-band metadata.
#[derive(Debug, Clone)]
pub struct BandInfo {
    pub name: String,
    pub data_type: String,
    pub nodata: Option<f64>,
}

/// Parsed raquet metadata (v0.4.0 format).
#[derive(Debug, Clone)]
pub struct RaquetMetadata {
    pub file_format: String,
    pub compression: String,
    pub compression_quality: i32,
    pub band_layout: String,
    pub block_width: i32,
    pub block_height: i32,
    pub min_zoom: i32,
    pub max_zoom: i32,
    pub pixel_zoom: i32,
    pub num_blocks: i32,
    pub scheme: String,
    pub crs: String,
    pub bands: Vec<(String, String)>,
    pub band_info: Vec<BandInfo>,
}

impl RaquetMetadata {
    pub fn is_interleaved(&self) -> bool {
        self.band_layout == "interleaved"
    }

    pub fn is_lossy_compression(&self) -> bool {
        self.compression == "jpeg" || self.compression == "webp"
    }

    pub fn num_bands(&self) -> usize {
        self.bands.len()
    }

    pub fn get_band_type(&self, band_index: usize) -> Option<&str> {
        self.bands.get(band_index).map(|(_, t)| t.as_str())
    }

    pub fn is_nodata(&self, band_index: usize, value: f64) -> bool {
        if let Some(info) = self.band_info.get(band_index) {
            if let Some(nd) = info.nodata {
                return value == nd || (value.is_nan() && nd.is_nan());
            }
        }
        false
    }
}

// ============================================================================
// Simple JSON extraction helpers
// ============================================================================

pub fn extract_json_string<'a>(json: &'a str, key: &str) -> &'a str {
    let search = format!("\"{}\":", key);
    let Some(pos) = json.find(&search) else {
        return "";
    };
    let rest = json[pos + search.len()..].trim_start_matches([' ', '\t', '\n']);
    if let Some(rest) = rest.strip_prefix('"') {
        let end = rest.find('"').unwrap_or(rest.len());
        &rest[..end]
    } else {
        // Numeric or other value – end at comma, } or ]
        let end = rest.find([',', '}', ']']).unwrap_or(rest.len());
        rest[..end].trim_end_matches([' ', '\t'])
    }
}

pub fn extract_json_int(json: &str, key: &str, default: i32) -> i32 {
    let v = extract_json_string(json, key);
    if v.is_empty() {
        return default;
    }
    v.parse().unwrap_or(default)
}

pub fn extract_json_double(json: &str, key: &str, default: f64) -> f64 {
    let v = extract_json_string(json, key);
    if v.is_empty() || v == "null" {
        return default;
    }
    v.parse().unwrap_or(default)
}

pub fn extract_json_has_value(json: &str, key: &str) -> bool {
    let v = extract_json_string(json, key);
    !v.is_empty() && v != "null"
}

/// Parse a nodata value, handling Zarr v3 string conventions.
pub fn parse_nodata_value(val: &str) -> f64 {
    match val {
        "" | "null" => 0.0,
        "NaN" => f64::NAN,
        "Infinity" => f64::INFINITY,
        "-Infinity" => f64::NEG_INFINITY,
        other => other.parse().unwrap_or(0.0),
    }
}

/// Extract a nested JSON object `"key": { ... }` as a string slice.
pub fn extract_json_object<'a>(json: &'a str, key: &str) -> &'a str {
    let search = format!("\"{}\":", key);
    let Some(pos) = json.find(&search) else {
        return "";
    };
    let rest = json[pos + search.len()..].trim_start_matches([' ', '\t', '\n']);
    if !rest.starts_with('{') {
        return "";
    }
    let mut depth = 1usize;
    let mut idx = 1;
    for c in rest[1..].chars() {
        match c {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    idx += 1;
                    break;
                }
            }
            _ => {}
        }
        idx += c.len_utf8();
    }
    &rest[..idx]
}

/// Parse the bands array from a metadata JSON string.
pub fn parse_bands_full(json: &str) -> (Vec<(String, String)>, Vec<BandInfo>) {
    let mut bands = Vec::new();
    let mut band_info = Vec::new();

    let Some(bands_pos) = json.find("\"bands\":") else {
        return (bands, band_info);
    };
    let rest = &json[bands_pos + 8..];
    let Some(arr_start) = rest.find('[') else {
        return (bands, band_info);
    };
    let rest = &rest[arr_start + 1..];
    let Some(arr_end) = rest.find(']') else {
        return (bands, band_info);
    };
    let bands_str = &rest[..arr_end];

    let mut pos = 0;
    while let Some(obj_start) = bands_str[pos..].find('{') {
        let obj_start = pos + obj_start;
        let Some(obj_end_rel) = bands_str[obj_start..].find('}') else {
            break;
        };
        let obj_end = obj_start + obj_end_rel + 1;
        let band_obj = &bands_str[obj_start..obj_end];

        let name = extract_json_string(band_obj, "name").to_string();
        let dtype = extract_json_string(band_obj, "type").to_string();
        if !name.is_empty() && !dtype.is_empty() {
            bands.push((name.clone(), dtype.clone()));
            let nodata = if extract_json_has_value(band_obj, "nodata") {
                let nd_str = extract_json_string(band_obj, "nodata");
                Some(parse_nodata_value(nd_str))
            } else {
                None
            };
            band_info.push(BandInfo {
                name,
                data_type: dtype,
                nodata,
            });
        }
        pos = obj_end;
    }
    (bands, band_info)
}

/// Parse a complete raquet metadata JSON string (v0.4.0 format).
pub fn parse_metadata(json: &str) -> RaquetMetadata {
    let file_format = extract_json_string(json, "file_format").to_string();
    let mut compression = extract_json_string(json, "compression").to_string();
    if compression.is_empty() {
        compression = "none".to_string();
    }
    let compression_quality = extract_json_int(json, "compression_quality", 0);
    let mut band_layout = extract_json_string(json, "band_layout").to_string();
    if band_layout.is_empty() {
        band_layout = "sequential".to_string();
    }
    let crs = extract_json_string(json, "crs").to_string();

    let tiling = extract_json_object(json, "tiling");
    let (min_zoom, max_zoom, pixel_zoom, num_blocks, block_width, block_height, scheme);
    if !tiling.is_empty() {
        min_zoom = extract_json_int(tiling, "min_zoom", 0);
        max_zoom = extract_json_int(tiling, "max_zoom", 26);
        pixel_zoom = extract_json_int(tiling, "pixel_zoom", 0);
        num_blocks = extract_json_int(tiling, "num_blocks", 0);
        block_width = extract_json_int(tiling, "block_width", 256);
        block_height = extract_json_int(tiling, "block_height", 256);
        scheme = extract_json_string(tiling, "scheme").to_string();
    } else {
        min_zoom = 0;
        max_zoom = 26;
        pixel_zoom = 0;
        num_blocks = 0;
        block_width = 256;
        block_height = 256;
        scheme = "quadbin".to_string();
    }

    let (bands, band_info) = parse_bands_full(json);

    RaquetMetadata {
        file_format,
        compression,
        compression_quality,
        band_layout,
        block_width,
        block_height,
        min_zoom,
        max_zoom,
        pixel_zoom,
        num_blocks,
        scheme,
        crs,
        bands,
        band_info,
    }
}
