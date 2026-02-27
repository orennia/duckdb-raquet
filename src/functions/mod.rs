pub mod macros;
/// DuckDB function registration modules.
pub mod quadbin_fns;
pub mod raster_fns;

use libduckdb_sys as ffi;
use std::ffi::CString;

/// Helper: register a scalar function set (with one or more overloads) on a connection.
///
/// `builders` is a slice of `(param_types, return_type, callback)`.
pub unsafe fn register_scalar_function_set(
    conn: ffi::duckdb_connection,
    name: &str,
    builders: &[(
        Vec<ffi::duckdb_logical_type>,
        ffi::duckdb_logical_type,
        ffi::duckdb_scalar_function_t,
    )],
) {
    let c_name = CString::new(name).unwrap();
    if builders.len() == 1 {
        let (params, ret, cb) = &builders[0];
        let func = ffi::duckdb_create_scalar_function();
        ffi::duckdb_scalar_function_set_name(func, c_name.as_ptr());
        for p in params {
            ffi::duckdb_scalar_function_add_parameter(func, *p);
        }
        ffi::duckdb_scalar_function_set_return_type(func, *ret);
        ffi::duckdb_scalar_function_set_function(func, *cb);
        ffi::duckdb_register_scalar_function(conn, func);
        let mut func = func;
        ffi::duckdb_destroy_scalar_function(&mut func);
    } else {
        let func_set = ffi::duckdb_create_scalar_function_set(c_name.as_ptr());
        for (params, ret, cb) in builders {
            let func = ffi::duckdb_create_scalar_function();
            for p in params {
                ffi::duckdb_scalar_function_add_parameter(func, *p);
            }
            ffi::duckdb_scalar_function_set_return_type(func, *ret);
            ffi::duckdb_scalar_function_set_function(func, *cb);
            ffi::duckdb_add_scalar_function_to_set(func_set, func);
            let mut func = func;
            ffi::duckdb_destroy_scalar_function(&mut func);
        }
        ffi::duckdb_register_scalar_function_set(conn, func_set);
        let mut func_set = func_set;
        ffi::duckdb_destroy_scalar_function_set(&mut func_set);
    }
}

/// Helper: create a DuckDB logical type for primitive types.
#[inline]
pub unsafe fn logical_type(t: ffi::DUCKDB_TYPE) -> ffi::duckdb_logical_type {
    ffi::duckdb_create_logical_type(t)
}

/// Helper: create a DuckDB list type wrapping a child type.
/// Consumes (destroys) the child type.
#[inline]
pub unsafe fn list_type(child: ffi::duckdb_logical_type) -> ffi::duckdb_logical_type {
    let lt = ffi::duckdb_create_list_type(child);
    let mut child = child;
    ffi::duckdb_destroy_logical_type(&mut child);
    lt
}

/// Helper: create a struct type from named members.
/// Destroys all member types after use.
pub unsafe fn struct_type(
    names: &[&str],
    types: &[ffi::duckdb_logical_type],
) -> ffi::duckdb_logical_type {
    let c_names: Vec<CString> = names.iter().map(|n| CString::new(*n).unwrap()).collect();
    let name_ptrs: Vec<*const std::os::raw::c_char> = c_names.iter().map(|n| n.as_ptr()).collect();
    let st = ffi::duckdb_create_struct_type(
        types.as_ptr() as *mut ffi::duckdb_logical_type,
        name_ptrs.as_ptr() as *mut *const std::os::raw::c_char,
        names.len() as ffi::idx_t,
    );
    for mut t in types.to_vec() {
        ffi::duckdb_destroy_logical_type(&mut t);
    }
    st
}

// ============================================================================
// Validity helpers
// ============================================================================

/// Check if row `i` is valid in a DuckDB validity bitmap.
/// A null bitmap pointer means all rows are valid.
#[inline]
pub unsafe fn row_is_valid(validity: *mut u64, i: ffi::idx_t) -> bool {
    if validity.is_null() {
        return true;
    }
    ffi::duckdb_validity_row_is_valid(validity, i)
}

/// Mark row `i` as NULL in the output vector.
#[inline]
pub unsafe fn set_row_invalid(output: ffi::duckdb_vector, i: ffi::idx_t) {
    ffi::duckdb_vector_ensure_validity_writable(output);
    let validity = ffi::duckdb_vector_get_validity(output);
    ffi::duckdb_validity_set_row_invalid(validity, i);
}

/// Read a null-terminated C string from a DuckDB varchar vector at row `i`.
///
/// # Safety
/// The caller must ensure `i` is in range and the row is valid.
pub unsafe fn read_varchar(vec: ffi::duckdb_vector, i: usize) -> String {
    let data = ffi::duckdb_vector_get_data(vec) as *const ffi::duckdb_string_t;
    let s = &*data.add(i);
    let len = s.value.inlined.length as usize;
    if len <= 12 {
        let raw = s.value.inlined.inlined.as_ptr() as *const u8;
        let bytes = std::slice::from_raw_parts(raw, len);
        std::str::from_utf8(bytes).unwrap_or("").to_string()
    } else {
        let ptr = s.value.pointer.ptr as *const u8;
        let bytes = std::slice::from_raw_parts(ptr, len);
        std::str::from_utf8(bytes).unwrap_or("").to_string()
    }
}

/// Read a BLOB (binary) value from a DuckDB vector at row `i`.
///
/// Returns a raw pointer and length. The caller is responsible for ensuring
/// the data is only accessed while the DuckDB vector remains valid (i.e.,
/// within the scalar function callback).
///
/// # Safety
/// The caller must ensure `i` is in range, the row is valid, and the returned
/// data is not used after the DuckDB vector is freed.
pub unsafe fn read_blob_raw(vec: ffi::duckdb_vector, i: usize) -> (*const u8, usize) {
    let data = ffi::duckdb_vector_get_data(vec) as *const ffi::duckdb_string_t;
    let s = &*data.add(i);
    let len = s.value.inlined.length as usize;
    if len <= 12 {
        let raw = s.value.inlined.inlined.as_ptr() as *const u8;
        (raw, len)
    } else {
        let ptr = s.value.pointer.ptr as *const u8;
        (ptr, len)
    }
}

/// Read a BLOB (binary) value from a DuckDB vector at row `i`.
///
/// # Safety
/// The caller must ensure `i` is in range, the row is valid, and the returned
/// slice is not used after the DuckDB vector is freed (i.e., only within the
/// scalar function callback that owns the vector).
pub unsafe fn read_blob<'a>(vec: ffi::duckdb_vector, i: usize) -> &'a [u8] {
    let (ptr, len) = read_blob_raw(vec, i);
    std::slice::from_raw_parts(ptr, len)
}

/// Write a string/blob result to a DuckDB output vector at row `i`.
pub unsafe fn write_string(output: ffi::duckdb_vector, i: usize, value: &[u8]) {
    ffi::duckdb_vector_assign_string_element_len(
        output,
        i as ffi::idx_t,
        value.as_ptr() as *const std::os::raw::c_char,
        value.len() as ffi::idx_t,
    );
}
