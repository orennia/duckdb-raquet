/// Band data decoding: gzip decompression and pixel value extraction.
/// Ported from the C++ band_decoder.cpp and band_decoder.hpp.
use flate2::read::GzDecoder;
use std::io::Read;

use crate::metadata::BandDataType;

// ============================================================================
// Gzip decompression
// ============================================================================

pub fn decompress_gzip(data: &[u8]) -> Result<Vec<u8>, String> {
    if data.is_empty() {
        return Ok(vec![]);
    }
    let mut decoder = GzDecoder::new(data);
    let mut result = Vec::new();
    decoder
        .read_to_end(&mut result)
        .map_err(|e| format!("gzip decompression failed: {}", e))?;
    Ok(result)
}

// ============================================================================
// float16 → f64 conversion
// ============================================================================

pub fn float16_to_f64(h: u16) -> f64 {
    let sign = (h >> 15) & 0x1;
    let exp = (h >> 10) & 0x1F;
    let mant = h & 0x3FF;
    if exp == 0 {
        if mant == 0 {
            return if sign != 0 { -0.0 } else { 0.0 };
        }
        let val = mant as f64 / 1024.0 * 2.0_f64.powi(-14);
        return if sign != 0 { -val } else { val };
    }
    if exp == 31 {
        return if mant == 0 {
            if sign != 0 {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            }
        } else {
            f64::NAN
        };
    }
    let val = (1.0 + mant as f64 / 1024.0) * 2.0_f64.powi(exp as i32 - 15);
    if sign != 0 {
        -val
    } else {
        val
    }
}

// ============================================================================
// Pixel extraction
// ============================================================================

/// Read a single pixel value from a decompressed pixel buffer.
pub fn get_pixel_value(data: &[u8], offset: usize, dtype: BandDataType) -> Result<f64, String> {
    let elem_size = dtype.byte_size();
    let byte_off = offset * elem_size;
    if byte_off + elem_size > data.len() {
        return Err(format!(
            "Pixel read at byte offset {} exceeds buffer of {} bytes",
            byte_off,
            data.len()
        ));
    }
    let val = match dtype {
        BandDataType::Uint8 => data[offset] as f64,
        BandDataType::Int8 => data[offset] as i8 as f64,
        BandDataType::Uint16 => {
            let v = u16::from_le_bytes(data[byte_off..byte_off + 2].try_into().unwrap());
            v as f64
        }
        BandDataType::Int16 => {
            let v = i16::from_le_bytes(data[byte_off..byte_off + 2].try_into().unwrap());
            v as f64
        }
        BandDataType::Uint32 => {
            let v = u32::from_le_bytes(data[byte_off..byte_off + 4].try_into().unwrap());
            v as f64
        }
        BandDataType::Int32 => {
            let v = i32::from_le_bytes(data[byte_off..byte_off + 4].try_into().unwrap());
            v as f64
        }
        BandDataType::Uint64 => {
            let v = u64::from_le_bytes(data[byte_off..byte_off + 8].try_into().unwrap());
            v as f64
        }
        BandDataType::Int64 => {
            let v = i64::from_le_bytes(data[byte_off..byte_off + 8].try_into().unwrap());
            v as f64
        }
        BandDataType::Float16 => {
            let v = u16::from_le_bytes(data[byte_off..byte_off + 2].try_into().unwrap());
            float16_to_f64(v)
        }
        BandDataType::Float32 => {
            let v = f32::from_le_bytes(data[byte_off..byte_off + 4].try_into().unwrap());
            v as f64
        }
        BandDataType::Float64 => {
            f64::from_le_bytes(data[byte_off..byte_off + 8].try_into().unwrap())
        }
    };
    Ok(val)
}

/// Decode a single pixel from band data.
pub fn decode_pixel(
    band_data: &[u8],
    dtype_str: &str,
    pixel_x: i32,
    pixel_y: i32,
    width: i32,
    compressed: bool,
) -> Result<f64, String> {
    let dtype = BandDataType::from_str(dtype_str)
        .ok_or_else(|| format!("Unknown data type: {}", dtype_str))?;
    let data = if compressed {
        decompress_gzip(band_data)?
    } else {
        band_data.to_vec()
    };
    let offset = (pixel_y * width + pixel_x) as usize;
    get_pixel_value(&data, offset, dtype)
}

/// Decode an entire band to a vector of f64 values.
pub fn decode_band(
    band_data: &[u8],
    dtype_str: &str,
    width: i32,
    height: i32,
    compressed: bool,
) -> Result<Vec<f64>, String> {
    let dtype = BandDataType::from_str(dtype_str)
        .ok_or_else(|| format!("Unknown data type: {}", dtype_str))?;
    let data = if compressed {
        decompress_gzip(band_data)?
    } else {
        band_data.to_vec()
    };
    let num_pixels = (width * height) as usize;
    let mut result = Vec::with_capacity(num_pixels);
    for i in 0..num_pixels {
        result.push(get_pixel_value(&data, i, dtype)?);
    }
    Ok(result)
}

/// Decode a pixel from interleaved (BIP) layout.
pub fn decode_pixel_interleaved(
    pixels_data: &[u8],
    dtype_str: &str,
    pixel_x: i32,
    pixel_y: i32,
    width: i32,
    band_index: usize,
    num_bands: usize,
    compression: &str,
) -> Result<f64, String> {
    let dtype = BandDataType::from_str(dtype_str)
        .ok_or_else(|| format!("Unknown data type: {}", dtype_str))?;
    let compressed = compression == "gzip";
    let data = if compressed {
        decompress_gzip(pixels_data)?
    } else {
        pixels_data.to_vec()
    };
    let pixel_linear = (pixel_y * width + pixel_x) as usize;
    let offset = pixel_linear * num_bands + band_index;
    get_pixel_value(&data, offset, dtype)
}

/// Decode an entire band from interleaved layout.
pub fn decode_band_interleaved(
    pixels_data: &[u8],
    dtype_str: &str,
    width: i32,
    height: i32,
    band_index: usize,
    num_bands: usize,
    compression: &str,
) -> Result<Vec<f64>, String> {
    let dtype = BandDataType::from_str(dtype_str)
        .ok_or_else(|| format!("Unknown data type: {}", dtype_str))?;
    let compressed = compression == "gzip";
    let data = if compressed {
        decompress_gzip(pixels_data)?
    } else {
        pixels_data.to_vec()
    };
    let num_pixels = (width * height) as usize;
    let mut result = Vec::with_capacity(num_pixels);
    for i in 0..num_pixels {
        let offset = i * num_bands + band_index;
        result.push(get_pixel_value(&data, offset, dtype)?);
    }
    Ok(result)
}

// ============================================================================
// Statistics
// ============================================================================

#[derive(Debug, Clone)]
pub struct BandStats {
    pub count: i64,
    pub sum: f64,
    pub mean: f64,
    pub min: f64,
    pub max: f64,
    pub stddev: f64,
}

/// Compute statistics for a band (streaming, no full pixel allocation).
pub fn compute_band_stats(
    band_data: &[u8],
    dtype_str: &str,
    width: i32,
    height: i32,
    compressed: bool,
    has_nodata: bool,
    nodata: f64,
) -> Result<BandStats, String> {
    let dtype = BandDataType::from_str(dtype_str)
        .ok_or_else(|| format!("Unknown data type: {}", dtype_str))?;
    let data = if compressed {
        decompress_gzip(band_data)?
    } else {
        band_data.to_vec()
    };
    let num_pixels = (width * height) as usize;

    let mut count: i64 = 0;
    let mut sum = 0.0_f64;
    let mut min = f64::MAX;
    let mut max = f64::MIN;
    let mut sum_sq = 0.0_f64;

    for i in 0..num_pixels {
        let v = get_pixel_value(&data, i, dtype)?;
        if has_nodata && (v == nodata || (v.is_nan() && nodata.is_nan())) {
            continue;
        }
        if v.is_nan() || v.is_infinite() {
            continue;
        }
        count += 1;
        sum += v;
        sum_sq += v * v;
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
    }

    let mean = if count > 0 { sum / count as f64 } else { 0.0 };
    let variance = if count > 0 {
        sum_sq / count as f64 - mean * mean
    } else {
        0.0
    };
    let stddev = variance.max(0.0).sqrt();
    let (min, max) = if count == 0 { (0.0, 0.0) } else { (min, max) };

    Ok(BandStats {
        count,
        sum,
        mean,
        min,
        max,
        stddev,
    })
}
