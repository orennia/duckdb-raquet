-- DuckDB Raquet Benchmark Queries
-- Run with: ./build/release/duckdb -c "LOAD raquet; .timer on; .read benchmark/queries_duckdb.sql"

-- ============================================================================
-- Configuration
-- ============================================================================
-- Test with GCS direct access (HTTPS URL)
-- For local files, replace with local path

.timer on

-- ============================================================================
-- Q1: Point Query - Get pixel value at a specific coordinate
-- ============================================================================
-- TCI.parquet - Sentinel-2 True Color
-- Location: Approx center of coverage area

SELECT '=== Q1: Point Query (ST_RasterValue) ===' as benchmark;

SELECT
    ST_RasterValue(block, band_1, ST_Point(33.5, 16.85), metadata) as red,
    ST_RasterValue(block, band_2, ST_Point(33.5, 16.85), metadata) as green,
    ST_RasterValue(block, band_3, ST_Point(33.5, 16.85), metadata) as blue
FROM read_raquet_at('https://storage.googleapis.com/sdsc_demo25/TCI.parquet', 33.5, 16.85);

-- ============================================================================
-- Q2: Tile Statistics - Compute stats for a single tile
-- ============================================================================

SELECT '=== Q2: Tile Statistics (ST_RasterSummaryStats) ===' as benchmark;

SELECT
    block,
    ST_RasterSummaryStats(band_1, metadata) as stats
FROM read_raquet('https://storage.googleapis.com/sdsc_demo25/TCI.parquet')
LIMIT 1;

-- ============================================================================
-- Q3: Region Statistics - Stats over a geographic polygon
-- ============================================================================

SELECT '=== Q3: Region Statistics ===' as benchmark;

SELECT
    count(*) as tile_count,
    sum((ST_RasterSummaryStats(band_1, metadata)).count) as total_pixels,
    avg((ST_RasterSummaryStats(band_1, metadata)).mean) as avg_mean
FROM read_raquet(
    'https://storage.googleapis.com/sdsc_demo25/TCI.parquet',
    'POLYGON((33.4 16.8, 33.6 16.8, 33.6 16.9, 33.4 16.9, 33.4 16.8))'::GEOMETRY
);

-- ============================================================================
-- Q4: Count tiles at each resolution
-- ============================================================================

SELECT '=== Q4: Resolution distribution ===' as benchmark;

SELECT
    quadbin_resolution(block) as resolution,
    count(*) as tile_count
FROM read_raquet('https://storage.googleapis.com/sdsc_demo25/TCI.parquet')
GROUP BY quadbin_resolution(block)
ORDER BY resolution;

-- ============================================================================
-- Q5: Bounding box query - Find all tiles in an area
-- ============================================================================

SELECT '=== Q5: Spatial filter ===' as benchmark;

SELECT count(*) as tiles_in_area
FROM read_raquet(
    'https://storage.googleapis.com/sdsc_demo25/TCI.parquet',
    'POLYGON((33.0 16.5, 34.0 16.5, 34.0 17.0, 33.0 17.0, 33.0 16.5))'::GEOMETRY
);

-- ============================================================================
-- Q6: Full table scan with aggregation
-- ============================================================================

SELECT '=== Q6: Full table aggregation ===' as benchmark;

SELECT
    count(*) as total_tiles,
    count(DISTINCT quadbin_resolution(block)) as resolution_levels,
    min(quadbin_resolution(block)) as min_resolution,
    max(quadbin_resolution(block)) as max_resolution
FROM read_raquet('https://storage.googleapis.com/sdsc_demo25/TCI.parquet');

SELECT '=== Benchmark Complete ===' as status;
