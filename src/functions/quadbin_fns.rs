//! Register all quadbin scalar functions with DuckDB.

use super::{
    list_type, logical_type, register_scalar_function_set, row_is_valid, set_row_invalid,
    struct_type,
};
use crate::quadbin::*;
use libduckdb_sys as ffi;

pub unsafe fn register_all(conn: ffi::duckdb_connection) {
    register_quadbin_from_tile(conn);
    register_quadbin_to_tile(conn);
    register_quadbin_from_lonlat(conn);
    register_quadbin_to_lonlat(conn);
    register_quadbin_resolution(conn);
    register_quadbin_to_bbox(conn);
    register_quadbin_to_bbox_mercator(conn);
    register_quadbin_to_parent(conn);
    register_quadbin_to_children(conn);
    register_quadbin_kring(conn);
    register_quadbin_siblings(conn);
    register_quadbin_polyfill(conn);
}

// ============================================================================
// quadbin_from_tile(x INT, y INT, z INT) -> UBIGINT
// ============================================================================

unsafe extern "C" fn quadbin_from_tile_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col_x = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_y = ffi::duckdb_data_chunk_get_vector(input, 1);
    let col_z = ffi::duckdb_data_chunk_get_vector(input, 2);
    let x_data = ffi::duckdb_vector_get_data(col_x) as *const i32;
    let y_data = ffi::duckdb_vector_get_data(col_y) as *const i32;
    let z_data = ffi::duckdb_vector_get_data(col_z) as *const i32;
    let vx = ffi::duckdb_vector_get_validity(col_x);
    let vy = ffi::duckdb_vector_get_validity(col_y);
    let vz = ffi::duckdb_vector_get_validity(col_z);
    let out_data = ffi::duckdb_vector_get_data(output) as *mut u64;
    for i in 0..n {
        if !row_is_valid(vx, i) || !row_is_valid(vy, i) || !row_is_valid(vz, i) {
            set_row_invalid(output, i);
            continue;
        }
        let x = *x_data.add(i as usize);
        let y = *y_data.add(i as usize);
        let z = *z_data.add(i as usize);
        match tile_to_cell(x, y, z) {
            Ok(cell) => *out_data.add(i as usize) = cell,
            Err(_) => set_row_invalid(output, i),
        }
    }
}

unsafe fn register_quadbin_from_tile(conn: ffi::duckdb_connection) {
    let int_t = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let int_t2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let int_t3 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let ret = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    register_scalar_function_set(
        conn,
        "quadbin_from_tile",
        &[(
            vec![int_t, int_t2, int_t3],
            ret,
            Some(quadbin_from_tile_impl),
        )],
    );
}

// ============================================================================
// quadbin_to_tile(cell UBIGINT) -> STRUCT(x INT, y INT, z INT)
// ============================================================================

unsafe extern "C" fn quadbin_to_tile_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let validity = ffi::duckdb_vector_get_validity(col);
    let x_child = ffi::duckdb_struct_vector_get_child(output, 0);
    let y_child = ffi::duckdb_struct_vector_get_child(output, 1);
    let z_child = ffi::duckdb_struct_vector_get_child(output, 2);
    let x_data = ffi::duckdb_vector_get_data(x_child) as *mut i32;
    let y_data = ffi::duckdb_vector_get_data(y_child) as *mut i32;
    let z_data = ffi::duckdb_vector_get_data(z_child) as *mut i32;
    for i in 0..n {
        if !row_is_valid(validity, i) {
            set_row_invalid(output, i);
            continue;
        }
        let cell = *in_data.add(i as usize);
        let (x, y, z) = cell_to_tile(cell);
        *x_data.add(i as usize) = x;
        *y_data.add(i as usize) = y;
        *z_data.add(i as usize) = z;
    }
}

unsafe fn register_quadbin_to_tile(conn: ffi::duckdb_connection) {
    let param = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let int_t = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let int_t2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let int_t3 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let ret = struct_type(&["x", "y", "z"], &[int_t, int_t2, int_t3]);
    register_scalar_function_set(
        conn,
        "quadbin_to_tile",
        &[(vec![param], ret, Some(quadbin_to_tile_impl))],
    );
}

// ============================================================================
// quadbin_from_lonlat(lon DOUBLE, lat DOUBLE, z INT) -> UBIGINT
// ============================================================================

unsafe extern "C" fn quadbin_from_lonlat_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col_lon = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_lat = ffi::duckdb_data_chunk_get_vector(input, 1);
    let col_z = ffi::duckdb_data_chunk_get_vector(input, 2);
    let lon_data = ffi::duckdb_vector_get_data(col_lon) as *const f64;
    let lat_data = ffi::duckdb_vector_get_data(col_lat) as *const f64;
    let z_data = ffi::duckdb_vector_get_data(col_z) as *const i32;
    let vlon = ffi::duckdb_vector_get_validity(col_lon);
    let vlat = ffi::duckdb_vector_get_validity(col_lat);
    let vz = ffi::duckdb_vector_get_validity(col_z);
    let out_data = ffi::duckdb_vector_get_data(output) as *mut u64;
    for i in 0..n {
        if !row_is_valid(vlon, i) || !row_is_valid(vlat, i) || !row_is_valid(vz, i) {
            set_row_invalid(output, i);
            continue;
        }
        let lon = *lon_data.add(i as usize);
        let lat = *lat_data.add(i as usize);
        let z = *z_data.add(i as usize);
        match lonlat_to_cell(lon, lat, z) {
            Ok(cell) => *out_data.add(i as usize) = cell,
            Err(_) => set_row_invalid(output, i),
        }
    }
}

unsafe fn register_quadbin_from_lonlat(conn: ffi::duckdb_connection) {
    let dbl = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let dbl2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let int_t = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let ret = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    register_scalar_function_set(
        conn,
        "quadbin_from_lonlat",
        &[(vec![dbl, dbl2, int_t], ret, Some(quadbin_from_lonlat_impl))],
    );
}

// ============================================================================
// quadbin_to_lonlat(cell UBIGINT) -> STRUCT(lon DOUBLE, lat DOUBLE)
// ============================================================================

unsafe extern "C" fn quadbin_to_lonlat_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let validity = ffi::duckdb_vector_get_validity(col);
    let lon_child = ffi::duckdb_struct_vector_get_child(output, 0);
    let lat_child = ffi::duckdb_struct_vector_get_child(output, 1);
    let lon_data = ffi::duckdb_vector_get_data(lon_child) as *mut f64;
    let lat_data = ffi::duckdb_vector_get_data(lat_child) as *mut f64;
    for i in 0..n {
        if !row_is_valid(validity, i) {
            set_row_invalid(output, i);
            continue;
        }
        let cell = *in_data.add(i as usize);
        let (lon, lat) = cell_to_lonlat(cell);
        *lon_data.add(i as usize) = lon;
        *lat_data.add(i as usize) = lat;
    }
}

unsafe fn register_quadbin_to_lonlat(conn: ffi::duckdb_connection) {
    let param = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let dbl = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let dbl2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let ret = struct_type(&["lon", "lat"], &[dbl, dbl2]);
    register_scalar_function_set(
        conn,
        "quadbin_to_lonlat",
        &[(vec![param], ret, Some(quadbin_to_lonlat_impl))],
    );
}

// ============================================================================
// quadbin_resolution(cell UBIGINT) -> INTEGER
// ============================================================================

unsafe extern "C" fn quadbin_resolution_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let validity = ffi::duckdb_vector_get_validity(col);
    let out_data = ffi::duckdb_vector_get_data(output) as *mut i32;
    for i in 0..n {
        if !row_is_valid(validity, i) {
            set_row_invalid(output, i);
            continue;
        }
        let cell = *in_data.add(i as usize);
        *out_data.add(i as usize) = cell_to_resolution(cell);
    }
}

unsafe fn register_quadbin_resolution(conn: ffi::duckdb_connection) {
    let param = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let ret = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    register_scalar_function_set(
        conn,
        "quadbin_resolution",
        &[(vec![param], ret, Some(quadbin_resolution_impl))],
    );
}

// ============================================================================
// quadbin_to_bbox(cell UBIGINT) -> STRUCT(xmin, ymin, xmax, ymax DOUBLE)
// ============================================================================

unsafe extern "C" fn quadbin_to_bbox_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let validity = ffi::duckdb_vector_get_validity(col);
    let c0 = ffi::duckdb_struct_vector_get_child(output, 0);
    let c1 = ffi::duckdb_struct_vector_get_child(output, 1);
    let c2 = ffi::duckdb_struct_vector_get_child(output, 2);
    let c3 = ffi::duckdb_struct_vector_get_child(output, 3);
    let d0 = ffi::duckdb_vector_get_data(c0) as *mut f64;
    let d1 = ffi::duckdb_vector_get_data(c1) as *mut f64;
    let d2 = ffi::duckdb_vector_get_data(c2) as *mut f64;
    let d3 = ffi::duckdb_vector_get_data(c3) as *mut f64;
    for i in 0..n {
        if !row_is_valid(validity, i) {
            set_row_invalid(output, i);
            continue;
        }
        let cell = *in_data.add(i as usize);
        let (tx, ty, tz) = cell_to_tile(cell);
        let (xmin, ymin, xmax, ymax) = tile_to_bbox_wgs84(tx, ty, tz);
        *d0.add(i as usize) = xmin;
        *d1.add(i as usize) = ymin;
        *d2.add(i as usize) = xmax;
        *d3.add(i as usize) = ymax;
    }
}

unsafe fn register_quadbin_to_bbox(conn: ffi::duckdb_connection) {
    let param = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let d0 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let d1 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let d2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let d3 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let ret = struct_type(&["xmin", "ymin", "xmax", "ymax"], &[d0, d1, d2, d3]);
    register_scalar_function_set(
        conn,
        "quadbin_to_bbox",
        &[(vec![param], ret, Some(quadbin_to_bbox_impl))],
    );
}

// ============================================================================
// quadbin_to_bbox_mercator(cell UBIGINT) -> STRUCT(xmin, ymin, xmax, ymax DOUBLE)
// ============================================================================

unsafe extern "C" fn quadbin_to_bbox_mercator_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let validity = ffi::duckdb_vector_get_validity(col);
    let c0 = ffi::duckdb_struct_vector_get_child(output, 0);
    let c1 = ffi::duckdb_struct_vector_get_child(output, 1);
    let c2 = ffi::duckdb_struct_vector_get_child(output, 2);
    let c3 = ffi::duckdb_struct_vector_get_child(output, 3);
    let d0 = ffi::duckdb_vector_get_data(c0) as *mut f64;
    let d1 = ffi::duckdb_vector_get_data(c1) as *mut f64;
    let d2 = ffi::duckdb_vector_get_data(c2) as *mut f64;
    let d3 = ffi::duckdb_vector_get_data(c3) as *mut f64;
    for i in 0..n {
        if !row_is_valid(validity, i) {
            set_row_invalid(output, i);
            continue;
        }
        let cell = *in_data.add(i as usize);
        let (tx, ty, tz) = cell_to_tile(cell);
        let (xmin, ymin, xmax, ymax) = tile_to_bbox_mercator(tx, ty, tz);
        *d0.add(i as usize) = xmin;
        *d1.add(i as usize) = ymin;
        *d2.add(i as usize) = xmax;
        *d3.add(i as usize) = ymax;
    }
}

unsafe fn register_quadbin_to_bbox_mercator(conn: ffi::duckdb_connection) {
    let param = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let d0 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let d1 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let d2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let d3 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_DOUBLE);
    let ret = struct_type(&["xmin", "ymin", "xmax", "ymax"], &[d0, d1, d2, d3]);
    register_scalar_function_set(
        conn,
        "quadbin_to_bbox_mercator",
        &[(vec![param], ret, Some(quadbin_to_bbox_mercator_impl))],
    );
}

// ============================================================================
// quadbin_to_parent overloads
// ============================================================================

unsafe extern "C" fn quadbin_to_parent_1_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let validity = ffi::duckdb_vector_get_validity(col);
    let out_data = ffi::duckdb_vector_get_data(output) as *mut u64;
    for i in 0..n {
        if !row_is_valid(validity, i) {
            set_row_invalid(output, i);
            continue;
        }
        let cell = *in_data.add(i as usize);
        match cell_to_parent(cell) {
            Ok(p) => *out_data.add(i as usize) = p,
            Err(_) => set_row_invalid(output, i),
        }
    }
}

unsafe extern "C" fn quadbin_to_parent_2_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_res = ffi::duckdb_data_chunk_get_vector(input, 1);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let res_data = ffi::duckdb_vector_get_data(col_res) as *const i32;
    let vc = ffi::duckdb_vector_get_validity(col);
    let vr = ffi::duckdb_vector_get_validity(col_res);
    let out_data = ffi::duckdb_vector_get_data(output) as *mut u64;
    for i in 0..n {
        if !row_is_valid(vc, i) || !row_is_valid(vr, i) {
            set_row_invalid(output, i);
            continue;
        }
        let cell = *in_data.add(i as usize);
        let res = *res_data.add(i as usize);
        match cell_to_parent_at(cell, res) {
            Ok(p) => *out_data.add(i as usize) = p,
            Err(_) => set_row_invalid(output, i),
        }
    }
}

unsafe fn register_quadbin_to_parent(conn: ffi::duckdb_connection) {
    let p1 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let r1 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let p2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let res2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let r2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    register_scalar_function_set(
        conn,
        "quadbin_to_parent",
        &[
            (vec![p1], r1, Some(quadbin_to_parent_1_impl)),
            (vec![p2, res2], r2, Some(quadbin_to_parent_2_impl)),
        ],
    );
}

// ============================================================================
// quadbin_to_children overloads
// ============================================================================

unsafe fn write_ubigint_list(
    output: ffi::duckdb_vector,
    i: usize,
    offset: &mut u64,
    items: &[u64],
) {
    let len = items.len() as u64;
    let list_entries = ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
    (*list_entries.add(i)).offset = *offset;
    (*list_entries.add(i)).length = len;

    ffi::duckdb_list_vector_reserve(output, *offset + len);
    let child = ffi::duckdb_list_vector_get_child(output);
    let child_data = ffi::duckdb_vector_get_data(child) as *mut u64;
    for (j, &v) in items.iter().enumerate() {
        *child_data.add((*offset + j as u64) as usize) = v;
    }
    *offset += len;
    ffi::duckdb_list_vector_set_size(output, *offset);
}

unsafe extern "C" fn quadbin_to_children_1_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let validity = ffi::duckdb_vector_get_validity(col);
    let mut offset = 0u64;
    for i in 0..n {
        if !row_is_valid(validity, i) {
            set_row_invalid(output, i);
            let list_entries = ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
            (*list_entries.add(i as usize)).offset = offset;
            (*list_entries.add(i as usize)).length = 0;
            continue;
        }
        let cell = *in_data.add(i as usize);
        match cell_to_children(cell) {
            Ok(children) => write_ubigint_list(output, i as usize, &mut offset, &children),
            Err(_) => {
                set_row_invalid(output, i);
                let list_entries =
                    ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
                (*list_entries.add(i as usize)).offset = offset;
                (*list_entries.add(i as usize)).length = 0;
            }
        }
    }
}

unsafe extern "C" fn quadbin_to_children_2_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_res = ffi::duckdb_data_chunk_get_vector(input, 1);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let res_data = ffi::duckdb_vector_get_data(col_res) as *const i32;
    let vc = ffi::duckdb_vector_get_validity(col);
    let vr = ffi::duckdb_vector_get_validity(col_res);
    let mut offset = 0u64;
    for i in 0..n {
        if !row_is_valid(vc, i) || !row_is_valid(vr, i) {
            set_row_invalid(output, i);
            let list_entries = ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
            (*list_entries.add(i as usize)).offset = offset;
            (*list_entries.add(i as usize)).length = 0;
            continue;
        }
        let cell = *in_data.add(i as usize);
        let res = *res_data.add(i as usize);
        match cell_to_children_at(cell, res) {
            Ok(children) => write_ubigint_list(output, i as usize, &mut offset, &children),
            Err(_) => {
                set_row_invalid(output, i);
                let list_entries =
                    ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
                (*list_entries.add(i as usize)).offset = offset;
                (*list_entries.add(i as usize)).length = 0;
            }
        }
    }
}

unsafe fn register_quadbin_to_children(conn: ffi::duckdb_connection) {
    let ubigint_t1 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let child_ret1 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let ret1 = list_type(child_ret1);

    let ubigint_t2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let int_t2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let child_ret2 = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let ret2 = list_type(child_ret2);

    register_scalar_function_set(
        conn,
        "quadbin_to_children",
        &[
            (vec![ubigint_t1], ret1, Some(quadbin_to_children_1_impl)),
            (
                vec![ubigint_t2, int_t2],
                ret2,
                Some(quadbin_to_children_2_impl),
            ),
        ],
    );
}

// ============================================================================
// quadbin_kring(cell UBIGINT, k INT) -> LIST(UBIGINT)
// ============================================================================

unsafe extern "C" fn quadbin_kring_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_k = ffi::duckdb_data_chunk_get_vector(input, 1);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let k_data = ffi::duckdb_vector_get_data(col_k) as *const i32;
    let vc = ffi::duckdb_vector_get_validity(col);
    let vk = ffi::duckdb_vector_get_validity(col_k);
    let mut offset = 0u64;
    for i in 0..n {
        if !row_is_valid(vc, i) || !row_is_valid(vk, i) {
            set_row_invalid(output, i);
            let list_entries = ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
            (*list_entries.add(i as usize)).offset = offset;
            (*list_entries.add(i as usize)).length = 0;
            continue;
        }
        let cell = *in_data.add(i as usize);
        let k = *k_data.add(i as usize);
        match cell_kring(cell, k) {
            Ok(neighbors) => write_ubigint_list(output, i as usize, &mut offset, &neighbors),
            Err(_) => {
                set_row_invalid(output, i);
                let list_entries =
                    ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
                (*list_entries.add(i as usize)).offset = offset;
                (*list_entries.add(i as usize)).length = 0;
            }
        }
    }
}

unsafe fn register_quadbin_kring(conn: ffi::duckdb_connection) {
    let p = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let k = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let child_ret = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let ret = list_type(child_ret);
    register_scalar_function_set(
        conn,
        "quadbin_kring",
        &[(vec![p, k], ret, Some(quadbin_kring_impl))],
    );
}

// ============================================================================
// quadbin_siblings(cell UBIGINT) -> LIST(UBIGINT)
// ============================================================================

unsafe extern "C" fn quadbin_siblings_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col = ffi::duckdb_data_chunk_get_vector(input, 0);
    let in_data = ffi::duckdb_vector_get_data(col) as *const u64;
    let validity = ffi::duckdb_vector_get_validity(col);
    let mut offset = 0u64;
    for i in 0..n {
        if !row_is_valid(validity, i) {
            set_row_invalid(output, i);
            let list_entries = ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
            (*list_entries.add(i as usize)).offset = offset;
            (*list_entries.add(i as usize)).length = 0;
            continue;
        }
        let cell = *in_data.add(i as usize);
        match cell_siblings(cell) {
            Ok(siblings) => write_ubigint_list(output, i as usize, &mut offset, &siblings),
            Err(_) => {
                set_row_invalid(output, i);
                let list_entries =
                    ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
                (*list_entries.add(i as usize)).offset = offset;
                (*list_entries.add(i as usize)).length = 0;
            }
        }
    }
}

unsafe fn register_quadbin_siblings(conn: ffi::duckdb_connection) {
    let p = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let child_ret = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let ret = list_type(child_ret);
    register_scalar_function_set(
        conn,
        "quadbin_siblings",
        &[(vec![p], ret, Some(quadbin_siblings_impl))],
    );
}

// ============================================================================
// QUADBIN_POLYFILL(geometry BLOB, resolution INT) -> LIST(UBIGINT)
// ============================================================================

unsafe extern "C" fn quadbin_polyfill_impl(
    _info: ffi::duckdb_function_info,
    input: ffi::duckdb_data_chunk,
    output: ffi::duckdb_vector,
) {
    let n = ffi::duckdb_data_chunk_get_size(input);
    let col_geom = ffi::duckdb_data_chunk_get_vector(input, 0);
    let col_res = ffi::duckdb_data_chunk_get_vector(input, 1);
    let vg = ffi::duckdb_vector_get_validity(col_geom);
    let vr = ffi::duckdb_vector_get_validity(col_res);
    let res_data = ffi::duckdb_vector_get_data(col_res) as *const i32;
    let mut offset = 0u64;
    for i in 0..n {
        if !row_is_valid(vg, i) || !row_is_valid(vr, i) {
            set_row_invalid(output, i);
            let list_entries = ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
            (*list_entries.add(i as usize)).offset = offset;
            (*list_entries.add(i as usize)).length = 0;
            continue;
        }
        let res = *res_data.add(i as usize);
        let geom = super::read_blob(col_geom, i as usize);
        match quadbin_polyfill(geom, res) {
            Ok(cells) => write_ubigint_list(output, i as usize, &mut offset, &cells),
            Err(_) => {
                set_row_invalid(output, i);
                let list_entries =
                    ffi::duckdb_vector_get_data(output) as *mut ffi::duckdb_list_entry;
                (*list_entries.add(i as usize)).offset = offset;
                (*list_entries.add(i as usize)).length = 0;
            }
        }
    }
}

unsafe fn register_quadbin_polyfill(conn: ffi::duckdb_connection) {
    let geom = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_BLOB);
    let res = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_INTEGER);
    let child_ret = logical_type(ffi::DUCKDB_TYPE_DUCKDB_TYPE_UBIGINT);
    let ret = list_type(child_ret);
    register_scalar_function_set(
        conn,
        "QUADBIN_POLYFILL",
        &[(vec![geom, res], ret, Some(quadbin_polyfill_impl))],
    );
}
