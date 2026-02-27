// Raquet DuckDB Extension – Rust port of the raster quadbin extension.
// Entry point follows the rusty-sheet pattern for C-API based extensions.

pub(crate) mod band_decoder;
pub(crate) mod functions;
pub(crate) mod metadata;
pub(crate) mod quadbin;

use duckdb::Connection;
use libduckdb_sys as ffi;

/// Minimum DuckDB version required for the C extension API.
/// Must match `TARGET_DUCKDB_VERSION` in the Makefile.
const DUCKDB_MIN_VERSION: &str = "v1.2.0";

/// Internal entry point (with error propagation).
///
/// # Safety
/// This function is called by DuckDB with valid extension info and access pointers.
pub unsafe fn raquet_init_c_api_internal(
    info: ffi::duckdb_extension_info,
    access: *const ffi::duckdb_extension_access,
) -> Result<bool, Box<dyn std::error::Error>> {
    // Initialize the DuckDB extension C API function pointers.
    let have_api_struct =
        unsafe { ffi::duckdb_rs_extension_api_init(info, access, DUCKDB_MIN_VERSION)? };
    if !have_api_struct {
        return Ok(false);
    }

    let db: ffi::duckdb_database = unsafe { *(*access).get_database.unwrap()(info) };

    // Open a high-level connection for registering SQL macros.
    let rust_conn = unsafe { Connection::open_from_raw(db.cast())? };

    // Open a second raw connection for the C-API scalar function registration.
    let mut raw_conn: ffi::duckdb_connection = std::ptr::null_mut();
    let state = unsafe { ffi::duckdb_connect(db, &mut raw_conn) };
    if state != ffi::duckdb_state_DuckDBSuccess {
        return Err("Failed to open raw DuckDB connection".into());
    }

    // Register all native scalar functions via the C API.
    unsafe {
        functions::quadbin_fns::register_all(raw_conn);
        functions::raster_fns::register_all(raw_conn);
    }

    // Register SQL macros via the high-level connection.
    functions::macros::register_all(&rust_conn)?;

    // Clean up the raw connection (functions are registered in the DB, not the connection).
    unsafe { ffi::duckdb_disconnect(&mut raw_conn) };

    Ok(true)
}

/// DuckDB extension C-API entry point.
///
/// DuckDB calls this function when loading an extension that exports
/// `{name}_init_c_api`.
///
/// # Safety
/// This function is called by DuckDB with valid extension info and access pointers.
#[unsafe(no_mangle)]
#[allow(clippy::not_unsafe_ptr_arg_deref)]
pub extern "C" fn raquet_init_c_api(
    info: ffi::duckdb_extension_info,
    access: *const ffi::duckdb_extension_access,
) -> bool {
    match unsafe { raquet_init_c_api_internal(info, access) } {
        Ok(result) => result,
        Err(e) => {
            let error_str = e.to_string();
            match std::ffi::CString::new(error_str) {
                Ok(c_err) => unsafe {
                    (*access).set_error.unwrap()(info, c_err.as_ptr());
                },
                Err(_) => unsafe {
                    let fallback = c"Raquet extension failed to initialise \
                        (and failed to allocate error string)";
                    (*access).set_error.unwrap()(info, fallback.as_ptr());
                },
            }
            false
        }
    }
}

// When loading by full path (e.g. LOAD 'target/debug/libraquet.so'), DuckDB
// uses the file's base name ("libraquet") to construct the function name.
// Export the alias so both `LOAD raquet` and `LOAD '.../libraquet.so'` work.
#[unsafe(no_mangle)]
pub extern "C" fn libraquet_init_c_api(
    info: ffi::duckdb_extension_info,
    access: *const ffi::duckdb_extension_access,
) -> bool {
    raquet_init_c_api(info, access)
}
