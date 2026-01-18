# DuckDB Raquet Extension

A DuckDB extension for working with [Raquet](https://github.com/CartoDB/raquet) raster data stored in Parquet format.

## Overview

Raquet is a specification for storing raster data in Apache Parquet with QUADBIN spatial indexing. This extension provides:

- **QUADBIN functions** for spatial indexing (tile ↔ cell ↔ lon/lat conversions)
- **Raster value extraction** (PostGIS-like `ST_RasterValue`)
- **Summary statistics** (PostGIS-like `ST_RasterSummaryStats`)
- **Band data decoding** for binary pixel data

## Installation

```sql
-- When published to community extensions:
INSTALL raquet FROM community;
LOAD raquet;
```

## Usage

### Reading Raquet Files

Raquet files are standard Parquet files, so you can read them directly:

```sql
-- Read from file
SELECT * FROM read_parquet('raster.parquet') WHERE block != 0 LIMIT 10;

-- Or create a table
CREATE TABLE dem AS SELECT * FROM read_parquet('elevation.parquet');
```

### QUADBIN Functions

```sql
-- Convert tile coordinates to QUADBIN cell
SELECT quadbin_from_tile(4096, 2048, 13);  -- Returns UBIGINT

-- Convert QUADBIN cell back to tile coordinates
SELECT (quadbin_to_tile(5234261499580514304)).*;  -- Returns {x, y, z}

-- Convert lon/lat to QUADBIN cell at resolution
SELECT quadbin_from_lonlat(-73.9857, 40.7484, 13);

-- Get center lon/lat of a QUADBIN cell
SELECT (quadbin_to_lonlat(5234261499580514304)).*;  -- Returns {lon, lat}

-- Get resolution from QUADBIN cell
SELECT quadbin_resolution(5234261499580514304);  -- Returns 13

-- Get bounding box of QUADBIN cell
SELECT (quadbin_to_bbox(5234261499580514304)).*;  -- Returns {min_lon, min_lat, max_lon, max_lat}

-- Get pixel coordinates within tile for a lon/lat
SELECT (quadbin_pixel_xy(-73.9857, 40.7484, 13, 256)).*;  -- Returns {pixel_x, pixel_y}
```

### Extracting Raster Values

```sql
-- Get pixel value at specific x,y within a tile
SELECT raquet_pixel(
    band_1,           -- BLOB: band data
    'int16',          -- VARCHAR: data type
    128,              -- INT: pixel x
    128,              -- INT: pixel y
    256,              -- INT: tile width
    'gzip'            -- VARCHAR: compression ('gzip' or 'none')
) FROM dem WHERE block = 5234261499580514304;

-- Get pixel value at lon/lat (PostGIS-style)
SELECT ST_RasterValue(
    block,            -- UBIGINT: quadbin cell
    band_1,           -- BLOB: band data
    -73.9857,         -- DOUBLE: longitude
    40.7484,          -- DOUBLE: latitude
    'int16',          -- VARCHAR: data type
    256,              -- INT: tile size
    'gzip'            -- VARCHAR: compression
) FROM dem WHERE block = quadbin_from_lonlat(-73.9857, 40.7484, 13);
```

### Joining Points with Raster

```sql
-- Get elevation at point locations
SELECT
    p.id,
    ST_RasterValue(r.block, r.band_1, ST_X(p.geom), ST_Y(p.geom), 'int16', 256, 'gzip') as elevation
FROM points p
JOIN dem r ON quadbin_from_lonlat(ST_X(p.geom), ST_Y(p.geom), 13) = r.block
WHERE r.block != 0;
```

### Summary Statistics

```sql
-- Get statistics for a tile
SELECT
    block,
    (ST_RasterSummaryStats(band_1, 'int16', 256, 256, 'gzip')).*
FROM dem
WHERE block != 0
LIMIT 5;

-- With nodata value filtering
SELECT
    block,
    (ST_RasterSummaryStats(band_1, 'int16', 256, 256, 'gzip', -9999)).*
FROM dem
WHERE block != 0;
```

### Decoding Full Bands

```sql
-- Decode entire band to array (useful for analysis)
SELECT
    block,
    raquet_decode_band(band_1, 'int16', 256, 256, 'gzip') as pixels
FROM dem
WHERE block != 0
LIMIT 1;
```

### Working with Metadata

```sql
-- Get the metadata row (block = 0)
SELECT metadata
FROM read_parquet('raster.parquet')
WHERE block = 0;

-- Helper functions
SELECT raquet_is_metadata_row(block), raquet_is_data_row(block)
FROM dem
LIMIT 5;
```

## Function Reference

### QUADBIN Functions

| Function | Description |
|----------|-------------|
| `quadbin_from_tile(x, y, z)` | Convert tile coordinates to QUADBIN cell |
| `quadbin_to_tile(cell)` | Convert QUADBIN cell to tile coordinates |
| `quadbin_from_lonlat(lon, lat, resolution)` | Convert lon/lat to QUADBIN cell |
| `quadbin_to_lonlat(cell)` | Get center lon/lat of QUADBIN cell |
| `quadbin_resolution(cell)` | Get resolution level from QUADBIN cell |
| `quadbin_to_bbox(cell)` | Get bounding box of QUADBIN cell |
| `quadbin_pixel_xy(lon, lat, resolution, tile_size)` | Get pixel coordinates within tile |

### Raster Functions

| Function | Description |
|----------|-------------|
| `raquet_pixel(band, dtype, x, y, width, compression)` | Get pixel value at x,y |
| `ST_RasterValue(block, band, lon, lat, dtype, width, compression)` | Get pixel value at lon/lat |
| `ST_RasterSummaryStats(band, dtype, width, height, compression[, nodata])` | Get tile statistics |
| `raquet_decode_band(band, dtype, width, height, compression)` | Decode band to array |

### Metadata Functions

| Function | Description |
|----------|-------------|
| `raquet_is_metadata_row(block)` | Check if row is metadata (block = 0) |
| `raquet_is_data_row(block)` | Check if row is data (block != 0) |

## Supported Data Types

- `uint8`, `int8`
- `uint16`, `int16`
- `uint32`, `int32`
- `uint64`, `int64`
- `float32`, `float64`

## Building from Source

```bash
# Clone the repository
git clone https://github.com/CartoDB/duckdb-raquet.git
cd duckdb-raquet

# Initialize submodules (DuckDB)
git submodule update --init --recursive

# Build
make release

# Run tests
make test
```

## Integration with DuckDB Spatial

This extension works seamlessly with the DuckDB spatial extension:

```sql
LOAD spatial;
LOAD raquet;

-- Create points table with geometry
CREATE TABLE points AS
SELECT 1 as id, ST_Point(-73.9857, 40.7484) as geom
UNION ALL
SELECT 2, ST_Point(-122.4194, 37.7749);

-- Join with raquet raster
SELECT
    p.id,
    ST_RasterValue(r.block, r.band_1, ST_X(p.geom), ST_Y(p.geom), 'int16', 256, 'gzip') as value
FROM points p
JOIN dem r ON quadbin_from_lonlat(ST_X(p.geom), ST_Y(p.geom), 13) = r.block
WHERE r.block != 0;
```

## License

Apache 2.0

## Related Projects

- [Raquet Specification](https://github.com/CartoDB/raquet) - The Raquet format specification
- [DuckDB Spatial](https://github.com/duckdb/duckdb-spatial) - Spatial extension for DuckDB
- [CARTO Analytics Toolbox](https://docs.carto.com/data-and-analysis/analytics-toolbox-for-bigquery/sql-reference/raster) - Raster functions for data warehouses
