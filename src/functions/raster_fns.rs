//! Register all raster scalar functions with DuckDB.

use super::{
    list_type, logical_type, read_blob, read_varchar, register_scalar_function_set, row_is_valid,
    set_row_invalid, struct_type, write_string,
};
use crate::band_decoder::*;
use crate::metadata::parse_metadata;
use crate::quadbin::{cell_to_tile, lonlat_to_pixel};
use libduckdb_sys as ffi;

pub unsafe fn register_all(conn: ffi::duckdb_connection) {
    register_raquet_pixel(conn);
    register_raquet_pixel_lonlat(conn);
    register_raquet_decode_band(conn);
    register_raquet_band_stats(conn);
    register_raquet_parse_metadata(conn);
    register_st_normalized_difference(conn);
}

// ============================================================================
// raquet_pixel(band BLOB, dtype VARCHAR, x INT, y INT, width INT, compression VARCHAR) -> DOUBLE
// ============================================================================

unsafe extern "C" fn raquet_pixel_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col_band = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_dtype = ffi::duckdb_data_chunk_get_vector(input, 1);
    let col_x = ffi::duckdb_data_chunk_get_vector(input, 2);
    let col_y = ffi::duckdb_data_chunk_get_vector(input, 3);
    let col_width = ffi::duckdb_data_chunk_get_vector(input, 4);
    let col_comp = ffi::duckdb_data_chunk_get_vector(input, 5);
    let vb = ffi::duckdb_vector_get_validity(col_band);
    let x_data = ffi::duckdb_vector_get_data(col_x) as *const i32;
    let y_data = ffi::duckdb_vector_get_data(col_y) as *const i32;
    let w_data = ffi::duckdb_vector_get_data(col_width) as *const i32;
    let out_data = ffi::duckdb_vector_get_data(output) as *mut f64;

    for i in 0..n {
        if !row_is_valid(vb, i) {
            set_row_invalid(output, i);
            continue;
        }
        let band = read_blob(col_band, i as usize);
        if band.is_empty() {
            set_row_invalid(output, i);
            continue;
        }
        let dtype = read_varchar(col_dtype, i as usize);
        let x = *x_data.add(i as usize);
        let y = *y_data.add(i as usize);
        let width = *w_data.add(i as usize);
        let comp = read_varchar(col_comp, i as usize);
        let compressed = comp == "gzip";
        if x < 0 || y < 0 || width <= 0 {
            set_row_invalid(output, i);
            continue;
        }
        match decode_pixel(band, &dtype, x, y, width, compressed) {
            Ok(v) => *out_data.add(i as usize) = v,
            Err(_) => set_row_invalid(output, i),
        }
    }
}

unsafe fn register_raquet_pixel(conn: ffi::duckdb_connection) {
    let params = vec![
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_BLOB),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR),
    ];
    let ret = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    register_scalar_function_set(
        conn,
        "raquet_pixel",
        &[(params.clone(), ret, Some(raquet_pixel_impl))],
    );
    let ret2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    register_scalar_function_set(conn, "ST_Value", &[(params, ret2, Some(raquet_pixel_impl))]);
}

// ============================================================================
// raquet_pixel_lonlat(band BLOB, dtype VARCHAR, lon DOUBLE, lat DOUBLE, cell UBIGINT, compression VARCHAR) -> DOUBLE
// ============================================================================

unsafe extern "C" fn raquet_pixel_lonlat_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col_band = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_dtype = ffi::duckdb_data_chunk_get_vector(input, 1);
    let col_lon = ffi::duckdb_data_chunk_get_vector(input, 2);
    let col_lat = ffi::duckdb_data_chunk_get_vector(input, 3);
    let col_cell = ffi::duckdb_data_chunk_get_vector(input, 4);
    let col_comp = ffi::duckdb_data_chunk_get_vector(input, 5);
    let vb = ffi::duckdb_vector_get_validity(col_band);
    let vc = ffi::duckdb_vector_get_validity(col_cell);
    let lon_data = ffi::duckdb_vector_get_data(col_lon) as *const f64;
    let lat_data = ffi::duckdb_vector_get_data(col_lat) as *const f64;
    let cell_data = ffi::duckdb_vector_get_data(col_cell) as *const u64;
    let out_data = ffi::duckdb_vector_get_data(output) as *mut f64;

    for i in 0..n {
        if !row_is_valid(vb, i) || !row_is_valid(vc, i) {
            set_row_invalid(output, i);
            continue;
        }
        let band = read_blob(col_band, i as usize);
        if band.is_empty() {
            set_row_invalid(output, i);
            continue;
        }
        let dtype = read_varchar(col_dtype, i as usize);
        let lon = *lon_data.add(i as usize);
        let lat = *lat_data.add(i as usize);
        let cell = *cell_data.add(i as usize);
        let comp = read_varchar(col_comp, i as usize);
        let compressed = comp == "gzip";

        let (_, _, z) = cell_to_tile(cell);
        let (px, py, _, _) = lonlat_to_pixel(lon, lat, z, 256);

        match decode_pixel(band, &dtype, px, py, 256, compressed) {
            Ok(v) => *out_data.add(i as usize) = v,
            Err(_) => set_row_invalid(output, i),
        }
    }
}

unsafe fn register_raquet_pixel_lonlat(conn: ffi::duckdb_connection) {
    let params = vec![
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_BLOB),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR),
    ];
    let ret = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    register_scalar_function_set(
        conn,
        "raquet_pixel_lonlat",
        &[(params, ret, Some(raquet_pixel_lonlat_impl))],
    );
}

// ============================================================================
// raquet_decode_band(band BLOB, dtype VARCHAR, width INT, height INT, compression VARCHAR) -> DOUBLE[]
// ============================================================================

unsafe extern "C" fn raquet_decode_band_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col_band = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_dtype = ffi::duckdb_data_chunk_get_vector(input, 1);
    let col_w = ffi::duckdb_data_chunk_get_vector(input, 2);
    let col_h = ffi::duckdb_data_chunk_get_vector(input, 3);
    let col_comp = ffi::duckdb_data_chunk_get_vector(input, 4);
    let vb = ffi::duckdb_vector_get_validity(col_band);
    let w_data = ffi::duckdb_vector_get_data(col_w) as *const i32;
    let h_data = ffi::duckdb_vector_get_data(col_h) as *const i32;
    let list_entries = ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
    let mut total_offset = 0u64;

    for i in 0..n {
        if !row_is_valid(vb, i) {
            (*list_entries.add(i as usize)).offset = total_offset;
            (*list_entries.add(i as usize)).length = 0;
            set_row_invalid(output, i);
            continue;
        }
        let band = read_blob(col_band, i as usize);
        if band.is_empty() {
            (*list_entries.add(i as usize)).offset = total_offset;
            (*list_entries.add(i as usize)).length = 0;
            set_row_invalid(output, i);
            continue;
        }
        let dtype = read_varchar(col_dtype, i as usize);
        let w = *w_data.add(i as usize);
        let h = *h_data.add(i as usize);
        let comp = read_varchar(col_comp, i as usize);
        let compressed = comp == "gzip";

        match decode_band(band, &dtype, w, h, compressed) {
            Ok(values) => {
                let len = values.len() as u64;
                (*list_entries.add(i as usize)).offset = total_offset;
                (*list_entries.add(i as usize)).length = len;
                ffi::duckdb_list_vector_reserve(output, total_offset + len);
                let child = ffi::duckdb_list_vector_get_child(output);
                let child_data = ffi::duckdb_vector_get_data(child) as *mut f64;
                for (j, &v) in values.iter().enumerate() {
                    *child_data.add((total_offset + j as u64) as usize) = v;
                }
                total_offset += len;
                ffi::duckdb_list_vector_set_size(output, total_offset);
            }
            Err(_) => {
                (*list_entries.add(i as usize)).offset = total_offset;
                (*list_entries.add(i as usize)).length = 0;
                set_row_invalid(output, i);
            }
        }
    }
}

unsafe fn register_raquet_decode_band(conn: ffi::duckdb_connection) {
    let params = vec![
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_BLOB),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR),
    ];
    let child_type = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let ret = list_type(child_type);
    register_scalar_function_set(
        conn,
        "raquet_decode_band",
        &[(params, ret, Some(raquet_decode_band_impl))],
    );
}

// ============================================================================
// raquet_band_stats / ST_BandStats / ST_RasterSummaryStats
// (band BLOB, dtype VARCHAR, width INT, height INT, compression VARCHAR, nodata DOUBLE)
// -> STRUCT(count BIGINT, sum DOUBLE, mean DOUBLE, min DOUBLE, max DOUBLE, stddev DOUBLE)
// ============================================================================

unsafe extern "C" fn raquet_band_stats_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col_band = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_dtype = ffi::duckdb_data_chunk_get_vector(input, 1);
    let col_w = ffi::duckdb_data_chunk_get_vector(input, 2);
    let col_h = ffi::duckdb_data_chunk_get_vector(input, 3);
    let col_comp = ffi::duckdb_data_chunk_get_vector(input, 4);
    let col_nodata = ffi::duckdb_data_chunk_get_vector(input, 5);
    let vb = ffi::duckdb_vector_get_validity(col_band);
    let vnd = ffi::duckdb_vector_get_validity(col_nodata);
    let w_data = ffi::duckdb_vector_get_data(col_w) as *const i32;
    let h_data = ffi::duckdb_vector_get_data(col_h) as *const i32;
    let nd_data = ffi::duckdb_vector_get_data(col_nodata) as *const f64;

    let c_count = ffi::duckdb_struct_vector_get_child(output, 0);
    let c_sum = ffi::duckdb_struct_vector_get_child(output, 1);
    let c_mean = ffi::duckdb_struct_vector_get_child(output, 2);
    let c_min = ffi::duckdb_struct_vector_get_child(output, 3);
    let c_max = ffi::duckdb_struct_vector_get_child(output, 4);
    let c_stddev = ffi::duckdb_struct_vector_get_child(output, 5);
    let d_count = ffi::duckdb_vector_get_data(c_count) as *mut i64;
    let d_sum = ffi::duckdb_vector_get_data(c_sum) as *mut f64;
    let d_mean = ffi::duckdb_vector_get_data(c_mean) as *mut f64;
    let d_min = ffi::duckdb_vector_get_data(c_min) as *mut f64;
    let d_max = ffi::duckdb_vector_get_data(c_max) as *mut f64;
    let d_stddev = ffi::duckdb_vector_get_data(c_stddev) as *mut f64;

    for i in 0..n {
        if !row_is_valid(vb, i) {
            set_row_invalid(output, i);
            continue;
        }
        let band = read_blob(col_band, i as usize);
        if band.is_empty() {
            set_row_invalid(output, i);
            continue;
        }
        let dtype = read_varchar(col_dtype, i as usize);
        let w = *w_data.add(i as usize);
        let h = *h_data.add(i as usize);
        let comp = read_varchar(col_comp, i as usize);
        let compressed = comp == "gzip";
        let has_nodata = row_is_valid(vnd, i);
        let nodata = if has_nodata {
            *nd_data.add(i as usize)
        } else {
            0.0
        };

        match compute_band_stats(band, &dtype, w, h, compressed, has_nodata, nodata) {
            Ok(stats) => {
                *d_count.add(i as usize) = stats.count;
                *d_sum.add(i as usize) = stats.sum;
                *d_mean.add(i as usize) = stats.mean;
                *d_min.add(i as usize) = stats.min;
                *d_max.add(i as usize) = stats.max;
                *d_stddev.add(i as usize) = stats.stddev;
            }
            Err(_) => set_row_invalid(output, i),
        }
    }
}

unsafe fn make_band_stats_types() -> (Vec<ffi::duckdb_logical_type>, ffi::duckdb_logical_type) {
    let params = vec![
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_BLOB),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE),
    ];
    let t_bigint = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_BIGINT);
    let t_dbl1 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let t_dbl2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let t_dbl3 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let t_dbl4 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let t_dbl5 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let ret = struct_type(
        &["count", "sum", "mean", "min", "max", "stddev"],
        &[t_bigint, t_dbl1, t_dbl2, t_dbl3, t_dbl4, t_dbl5],
    );
    (params, ret)
}

unsafe fn register_raquet_band_stats(conn: ffi::duckdb_connection) {
    let (params, ret) = make_band_stats_types();
    register_scalar_function_set(
        conn,
        "raquet_band_stats",
        &[(params, ret, Some(raquet_band_stats_impl))],
    );
    let (params2, ret2) = make_band_stats_types();
    register_scalar_function_set(
        conn,
        "ST_BandStats",
        &[(params2, ret2, Some(raquet_band_stats_impl))],
    );
    let (params3, ret3) = make_band_stats_types();
    register_scalar_function_set(
        conn,
        "ST_RasterSummaryStats",
        &[(params3, ret3, Some(raquet_band_stats_impl))],
    );
}

// ============================================================================
// raquet_parse_metadata(metadata_json VARCHAR)
// -> STRUCT(compression, compression_quality, band_layout, block_width, block_height,
//            min_zoom, max_zoom, num_bands)
// ============================================================================

unsafe extern "C" fn raquet_parse_metadata_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let validity = ffi::duckdb_vector_get_validity(col);

    let c_comp = ffi::duckdb_struct_vector_get_child(output, 0);
    let c_comp_q = ffi::duckdb_struct_vector_get_child(output, 1);
    let c_band_layout = ffi::duckdb_struct_vector_get_child(output, 2);
    let c_bw = ffi::duckdb_struct_vector_get_child(output, 3);
    let c_bh = ffi::duckdb_struct_vector_get_child(output, 4);
    let c_min_zoom = ffi::duckdb_struct_vector_get_child(output, 5);
    let c_max_zoom = ffi::duckdb_struct_vector_get_child(output, 6);
    let c_num_bands = ffi::duckdb_struct_vector_get_child(output, 7);

    let d_comp_q = ffi::duckdb_vector_get_data(c_comp_q) as *mut i32;
    let d_bw = ffi::duckdb_vector_get_data(c_bw) as *mut i32;
    let d_bh = ffi::duckdb_vector_get_data(c_bh) as *mut i32;
    let d_min_zoom = ffi::duckdb_vector_get_data(c_min_zoom) as *mut i32;
    let d_max_zoom = ffi::duckdb_vector_get_data(c_max_zoom) as *mut i32;
    let d_num_bands = ffi::duckdb_vector_get_data(c_num_bands) as *mut i32;

    for i in 0..n {
        if !row_is_valid(validity, i) {
            set_row_invalid(output, i);
            continue;
        }
        let json = read_varchar(col, i as usize);
        if json.is_empty() {
            set_row_invalid(output, i);
            continue;
        }
        let meta = parse_metadata(&json);
        write_string(c_comp, i as usize, meta.compression.as_bytes());
        *d_comp_q.add(i as usize) = meta.compression_quality;
        write_string(c_band_layout, i as usize, meta.band_layout.as_bytes());
        *d_bw.add(i as usize) = meta.block_width;
        *d_bh.add(i as usize) = meta.block_height;
        *d_min_zoom.add(i as usize) = meta.min_zoom;
        *d_max_zoom.add(i as usize) = meta.max_zoom;
        *d_num_bands.add(i as usize) = meta.num_bands() as i32;
    }
}

unsafe fn register_raquet_parse_metadata(conn: ffi::duckdb_connection) {
    let param = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR);
    let t_varchar1 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR);
    let t_int1 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let t_varchar2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR);
    let t_int2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let t_int3 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let t_int4 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let t_int5 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let t_int6 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let ret = struct_type(
        &[
            "compression",
            "compression_quality",
            "band_layout",
            "block_width",
            "block_height",
            "min_zoom",
            "max_zoom",
            "num_bands",
        ],
        &[
            t_varchar1, t_int1, t_varchar2, t_int2, t_int3, t_int4, t_int5, t_int6,
        ],
    );
    register_scalar_function_set(
        conn,
        "raquet_parse_metadata",
        &[(vec![param], ret, Some(raquet_parse_metadata_impl))],
    );
}

// ============================================================================
// ST_NormalizedDifference(band1 BLOB, band2 BLOB, metadata VARCHAR) -> DOUBLE[]
// ============================================================================

unsafe extern "C" fn st_normalized_difference_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col_b1 = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_b2 = ffi::duckdb_data_chunk_get_vector(input, 1);
    let col_meta = ffi::duckdb_data_chunk_get_vector(input, 2);
    let vb1 = ffi::duckdb_vector_get_validity(col_b1);
    let vb2 = ffi::duckdb_vector_get_validity(col_b2);
    let list_entries = ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
    let mut total_offset = 0u64;

    for i in 0..n {
        if !row_is_valid(vb1, i) || !row_is_valid(vb2, i) {
            (*list_entries.add(i as usize)).offset = total_offset;
            (*list_entries.add(i as usize)).length = 0;
            set_row_invalid(output, i);
            continue;
        }
        let b1 = read_blob(col_b1, i as usize);
        let b2 = read_blob(col_b2, i as usize);
        let meta_str = read_varchar(col_meta, i as usize);
        let meta = parse_metadata(&meta_str);
        let dtype = meta
            .bands
            .first()
            .map(|(_, t)| t.as_str())
            .unwrap_or("float32");
        let compressed = meta.compression == "gzip";
        let width = meta.block_width;
        let height = meta.block_height;

        let vals1 = match decode_band(b1, dtype, width, height, compressed) {
            Ok(v) => v,
            Err(_) => {
                (*list_entries.add(i as usize)).offset = total_offset;
                (*list_entries.add(i as usize)).length = 0;
                set_row_invalid(output, i);
                continue;
            }
        };
        let vals2 = match decode_band(b2, dtype, width, height, compressed) {
            Ok(v) => v,
            Err(_) => {
                (*list_entries.add(i as usize)).offset = total_offset;
                (*list_entries.add(i as usize)).length = 0;
                set_row_invalid(output, i);
                continue;
            }
        };

        let len = vals1.len().min(vals2.len()) as u64;
        (*list_entries.add(i as usize)).offset = total_offset;
        (*list_entries.add(i as usize)).length = len;
        ffi::duckdb_list_vector_reserve(output, total_offset + len);
        let child = ffi::duckdb_list_vector_get_child(output);
        let child_data = ffi::duckdb_vector_get_data(child) as *mut f64;
        for j in 0..len as usize {
            let v1 = vals1[j];
            let v2 = vals2[j];
            let nd = if (v1 + v2).abs() < 1e-15 {
                f64::NAN
            } else {
                (v1 - v2) / (v1 + v2)
            };
            *child_data.add((total_offset + j as u64) as usize) = nd;
        }
        total_offset += len;
        ffi::duckdb_list_vector_set_size(output, total_offset);
    }
}

unsafe fn register_st_normalized_difference(conn: ffi::duckdb_connection) {
    let params = vec![
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_BLOB),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_BLOB),
        logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_VARCHAR),
    ];
    let child_type = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let ret = list_type(child_type);
    register_scalar_function_set(
        conn,
        "ST_NormalizedDifference",
        &[(params, ret, Some(st_normalized_difference_impl))],
    );
}
