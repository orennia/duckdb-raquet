# Agent Instructions for duckdb-raquet

## Project Overview

`duckdb-raquet` is a C++ DuckDB extension for working with [Raquet](https://github.com/CartoDB/raquet) raster data stored in Apache Parquet format with QUADBIN spatial indexing. It exposes SQL functions for spatial indexing, pixel extraction, raster analytics, and cloud-native querying of remote Parquet files.

## Repository Structure

```
src/
  include/               # Header files
  quadbin/               # QUADBIN spatial index functions
  raster/                # Raster processing (band decoding, stats, pixel ops)
  metadata/              # Raquet metadata parsing
  table_functions/       # DuckDB table macro registrations
  raquet_extension.cpp   # Extension entry point & macro definitions
test/sql/                # SQL-based integration tests (DuckDB test framework)
.github/workflows/       # CI pipelines
CMakeLists.txt           # CMake build configuration
Makefile                 # Developer convenience targets
extension_config.cmake   # DuckDB extension build config
```

## Build System

The project uses CMake via DuckDB's extension build system.

**Prerequisites:**
- CMake 3.12+
- C++17 compatible compiler
- zlib (required)
- libjpeg (optional, enables JPEG lossy compression)
- libwebp (optional, enables WebP lossy compression)
- DuckDB source tree cloned as `./duckdb/` (run `make setup` to clone it)

**Common build commands:**
```bash
make setup      # Clone DuckDB source if not present
make release    # Build release binary
make debug      # Build debug binary
make test_sql   # Run SQL integration tests
make clean      # Remove build artifacts
make format     # Auto-format C++ source with clang-format
```

## Testing

Tests live in `test/sql/` and use the DuckDB SQL test framework (`.test` files). Run them with:

```bash
make test_sql
```

Each `.test` file contains SQL statements and expected outputs. To add a test, create or edit a `.test` file in `test/sql/`. Follow the naming convention of existing files (e.g., `quadbin.test`, `raster.test`).

## Code Conventions

- **Language:** C++17
- **Namespace:** All code lives in the `duckdb` namespace
- **Formatting:** `clang-format` via `make format`; match the style of existing files
- **Function registration:** New SQL functions are registered in `raquet_extension.cpp` via `Register*Functions()` helpers defined per module
- **Headers:** Declare new public functions in `src/include/`
- **No magic numbers:** Use named constants for tile sizes, resolutions, etc.
- **Error handling:** Use `duckdb::InvalidInputException` or `duckdb::NotImplementedException` for user-facing errors

## Key Concepts

- **QUADBIN** – 64-bit integer encoding of Web Mercator tile coordinates (x, y, zoom). All spatial operations reference these identifiers.
- **Raquet v0.3.0+** – Parquet files where `block=0` row holds metadata JSON; data rows have `block != 0`.
- **Band data** – Pixel values stored as compressed BLOBs (`gzip` or uncompressed by default; `JPEG`/`WebP` with v0.4.0).
- **read_raquet** – Table macro that wraps `read_parquet`, propagates metadata, and supports optional spatial filtering.

## Adding New SQL Functions

1. Implement the function in the appropriate module under `src/raster/`, `src/quadbin/`, etc.
2. Register it via a `Register*Functions(ExtensionLoader &loader)` call in `raquet_extension.cpp`.
3. Add or update `.test` files in `test/sql/` covering the new behavior.
4. Export any new public headers in `src/include/`.

## CI

CI runs on push/PR to `main`/`master` via `.github/workflows/ci.yml`:
- Builds on Linux (ubuntu-latest) and macOS (macos-latest)
- Runs SQL tests with `make test_sql`
- Uploads the compiled `.duckdb_extension` as a build artifact
