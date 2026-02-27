.PHONY: all clean debug release test configure

PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
EXTENSION_NAME := raquet

# Set to 1 to enable Unstable API (binaries will only work on TARGET_DUCKDB_VERSION)
USE_UNSTABLE_C_API := 0

# Target DuckDB version
TARGET_DUCKDB_VERSION := v1.2.0

all: debug

# Configure: ensure Rust toolchain is present
configure:
@which cargo > /dev/null || (echo "Rust toolchain not found. Install from https://rustup.rs/" && exit 1)
@echo "Rust toolchain OK: $$(cargo --version)"

# Debug build
debug: configure
cargo build

# Release build
release: configure
cargo build --release

# Run unit tests (Rust)
test: configure
cargo test

# Run SQL tests (requires a built extension and duckdb binary)
test_sql: debug
@echo "SQL tests require a DuckDB binary with the extension loaded."
@echo "Build the extension with 'make debug' and load it manually."

# Clean
clean:
cargo clean

# Format code
format:
cargo fmt

# Lint
lint:
cargo clippy
