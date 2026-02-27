/// Register SQL table macros with DuckDB via connection.execute_batch.
use duckdb::Connection;

pub fn register_all(conn: &Connection) -> Result<(), duckdb::Error> {
    // ========================================================================
    // read_raquet - read a raquet parquet file, propagating metadata
    // ========================================================================

    // 1-arg: basic read
    conn.execute_batch(
        "CREATE OR REPLACE MACRO read_raquet(file) AS TABLE
         WITH src AS (SELECT * FROM read_parquet(file))
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block != 0;",
    )?;

    // 2-arg: spatial filter with auto-detected resolution
    conn.execute_batch(
        "CREATE OR REPLACE MACRO read_raquet(file, geometry) AS TABLE
         WITH src AS (SELECT * FROM read_parquet(file)),
         file_resolution AS (
             SELECT quadbin_resolution(block) AS res
             FROM src
             WHERE block != 0
             LIMIT 1
         )
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block IN (
             SELECT UNNEST(QUADBIN_POLYFILL(geometry, (SELECT res FROM file_resolution)))
         )
         AND block != 0;",
    )?;

    // 3-arg: spatial filter with explicit resolution
    conn.execute_batch(
        "CREATE OR REPLACE MACRO read_raquet(file, geometry, resolution) AS TABLE
         WITH src AS (SELECT * FROM read_parquet(file))
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block IN (SELECT UNNEST(QUADBIN_POLYFILL(geometry, resolution)))
           AND block != 0;",
    )?;

    // ========================================================================
    // read_raquet_metadata - return only the metadata row
    // ========================================================================
    conn.execute_batch(
        "CREATE OR REPLACE MACRO read_raquet_metadata(file) AS TABLE
         SELECT metadata
         FROM read_parquet(file)
         WHERE block = 0;",
    )?;

    // ========================================================================
    // read_raquet_at - point query convenience function
    // ========================================================================

    // 2-arg: point geometry with auto-detected resolution
    conn.execute_batch(
        "CREATE OR REPLACE MACRO read_raquet_at(file, point) AS TABLE
         WITH src AS (SELECT * FROM read_parquet(file)),
         file_resolution AS (
             SELECT quadbin_resolution(block) AS res
             FROM src
             WHERE block != 0
             LIMIT 1
         )
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block = quadbin_from_lonlat(ST_X(point), ST_Y(point), (SELECT res FROM file_resolution));",
    )?;

    // 3-arg: lon/lat with auto-detected resolution
    conn.execute_batch(
        "CREATE OR REPLACE MACRO read_raquet_at(file, lon, lat) AS TABLE
         WITH src AS (SELECT * FROM read_parquet(file)),
         file_resolution AS (
             SELECT quadbin_resolution(block) AS res
             FROM src
             WHERE block != 0
             LIMIT 1
         )
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block = quadbin_from_lonlat(lon, lat, (SELECT res FROM file_resolution));",
    )?;

    // 4-arg: lon/lat with explicit resolution
    conn.execute_batch(
        "CREATE OR REPLACE MACRO read_raquet_at(file, lon, lat, resolution) AS TABLE
         WITH src AS (SELECT * FROM read_parquet(file))
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block = quadbin_from_lonlat(lon, lat, resolution);",
    )?;

    // ========================================================================
    // ST_Raster - read from a table (iceberg or any table)
    // ========================================================================

    // 1-arg: basic read from table
    conn.execute_batch(
        "CREATE OR REPLACE MACRO ST_Raster(tbl) AS TABLE
         WITH src AS (SELECT * FROM query_table(tbl))
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block != 0;",
    )?;

    // 2-arg: spatial filter with auto-detected resolution
    conn.execute_batch(
        "CREATE OR REPLACE MACRO ST_Raster(tbl, geometry) AS TABLE
         WITH src AS (SELECT * FROM query_table(tbl)),
         table_resolution AS (
             SELECT quadbin_resolution(block::UBIGINT) AS res
             FROM src
             WHERE block != 0
             LIMIT 1
         )
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block::UBIGINT IN (
             SELECT UNNEST(QUADBIN_POLYFILL(geometry, (SELECT res FROM table_resolution)))
         )
         AND block != 0;",
    )?;

    // 3-arg: spatial filter with explicit resolution
    conn.execute_batch(
        "CREATE OR REPLACE MACRO ST_Raster(tbl, geometry, resolution) AS TABLE
         WITH src AS (SELECT * FROM query_table(tbl))
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block::UBIGINT IN (SELECT UNNEST(QUADBIN_POLYFILL(geometry, resolution)))
           AND block != 0;",
    )?;

    // ========================================================================
    // ST_RasterAt - point query on a table
    // ========================================================================

    // 2-arg: point geometry with auto-detected resolution
    conn.execute_batch(
        "CREATE OR REPLACE MACRO ST_RasterAt(tbl, point) AS TABLE
         WITH src AS (SELECT * FROM query_table(tbl)),
         table_resolution AS (
             SELECT quadbin_resolution(block::UBIGINT) AS res
             FROM src
             WHERE block != 0
             LIMIT 1
         )
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block::UBIGINT = quadbin_from_lonlat(ST_X(point), ST_Y(point), (SELECT res FROM table_resolution));",
    )?;

    // 3-arg: lon/lat with auto-detected resolution
    conn.execute_batch(
        "CREATE OR REPLACE MACRO ST_RasterAt(tbl, lon, lat) AS TABLE
         WITH src AS (SELECT * FROM query_table(tbl)),
         table_resolution AS (
             SELECT quadbin_resolution(block::UBIGINT) AS res
             FROM src
             WHERE block != 0
             LIMIT 1
         )
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block::UBIGINT = quadbin_from_lonlat(lon, lat, (SELECT res FROM table_resolution));",
    )?;

    // 4-arg: lon/lat with explicit resolution
    conn.execute_batch(
        "CREATE OR REPLACE MACRO ST_RasterAt(tbl, lon, lat, resolution) AS TABLE
         WITH src AS (SELECT * FROM query_table(tbl))
         SELECT * REPLACE (
             (SELECT metadata FROM src WHERE block = 0) AS metadata
         )
         FROM src
         WHERE block::UBIGINT = quadbin_from_lonlat(lon, lat, resolution);",
    )?;

    Ok(())
}
