#[no_mangle]
pub extern "C" fn raquet_rust_version() -> *const std::ffi::c_char {
    c"0.1.0".as_ptr()
}
