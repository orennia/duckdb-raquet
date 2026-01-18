#pragma once

#include <cstdint>
#include <cmath>
#include <stdexcept>

namespace duckdb {
namespace quadbin {

// Constants
constexpr int MAX_RESOLUTION = 26;
constexpr double EARTH_RADIUS = 6378137.0;
constexpr double MAX_LATITUDE = 85.051128779806604;
constexpr double PI = 3.14159265358979323846;

// QUADBIN encoding constants (from quadbin Python library)
constexpr uint64_t HEADER = 0x4000000000000000ULL;
constexpr uint64_t MODE = 0x0800000000000000ULL;  // Mode 1 for cells
constexpr uint64_t FOOTER = 0x000FFFFFFFFFFFFFULL;

// Magic numbers for Morton code interleaving
constexpr uint64_t B0 = 0x5555555555555555ULL;
constexpr uint64_t B1 = 0x3333333333333333ULL;
constexpr uint64_t B2 = 0x0F0F0F0F0F0F0F0FULL;
constexpr uint64_t B3 = 0x00FF00FF00FF00FFULL;
constexpr uint64_t B4 = 0x0000FFFF0000FFFFULL;
constexpr uint64_t B5 = 0x00000000FFFFFFFFULL;

// Extract resolution from quadbin cell
inline int cell_to_resolution(uint64_t cell) {
    return static_cast<int>((cell >> 52) & 0x1F);
}

// Convert tile coordinates to quadbin cell
inline uint64_t tile_to_cell(int x, int y, int z) {
    if (z < 0 || z > MAX_RESOLUTION) {
        throw std::invalid_argument("Resolution must be between 0 and 26");
    }

    // Shift x and y to 32-bit range based on resolution
    uint64_t ux = static_cast<uint64_t>(x) << (32 - z);
    uint64_t uy = static_cast<uint64_t>(y) << (32 - z);

    // Morton code interleaving (spread bits)
    ux = (ux | (ux << 16)) & B4;
    ux = (ux | (ux << 8)) & B3;
    ux = (ux | (ux << 4)) & B2;
    ux = (ux | (ux << 2)) & B1;
    ux = (ux | (ux << 1)) & B0;

    uy = (uy | (uy << 16)) & B4;
    uy = (uy | (uy << 8)) & B3;
    uy = (uy | (uy << 4)) & B2;
    uy = (uy | (uy << 2)) & B1;
    uy = (uy | (uy << 1)) & B0;

    // Combine: HEADER | MODE | resolution | interleaved_xy | trailing_ones
    return HEADER | MODE | (static_cast<uint64_t>(z) << 52) |
           ((ux | (uy << 1)) >> 12) | (FOOTER >> (z * 2));
}

// Convert quadbin cell to tile coordinates
inline void cell_to_tile(uint64_t cell, int &x, int &y, int &z) {
    z = cell_to_resolution(cell);

    // Extract the quadkey portion and shift to working position
    uint64_t q = (cell & FOOTER) << 12;

    // De-interleave x (even bits) and y (odd bits)
    uint64_t ux = q;
    uint64_t uy = q >> 1;

    // Compact the bits back (reverse Morton code)
    ux = ux & B0;
    uy = uy & B0;

    ux = (ux | (ux >> 1)) & B1;
    uy = (uy | (uy >> 1)) & B1;

    ux = (ux | (ux >> 2)) & B2;
    uy = (uy | (uy >> 2)) & B2;

    ux = (ux | (ux >> 4)) & B3;
    uy = (uy | (uy >> 4)) & B3;

    ux = (ux | (ux >> 8)) & B4;
    uy = (uy | (uy >> 8)) & B4;

    ux = (ux | (ux >> 16)) & B5;
    uy = (uy | (uy >> 16)) & B5;

    // Shift back based on resolution
    x = static_cast<int>(ux >> (32 - z));
    y = static_cast<int>(uy >> (32 - z));
}

// Convert longitude/latitude to tile coordinates at given resolution
inline void lonlat_to_tile(double lon, double lat, int z, int &x, int &y) {
    // Clamp latitude to valid range
    if (lat > MAX_LATITUDE) lat = MAX_LATITUDE;
    if (lat < -MAX_LATITUDE) lat = -MAX_LATITUDE;

    // Convert to tile coordinates
    double n = std::pow(2.0, z);
    x = static_cast<int>(std::floor((lon + 180.0) / 360.0 * n));

    double lat_rad = lat * PI / 180.0;
    y = static_cast<int>(std::floor((1.0 - std::asinh(std::tan(lat_rad)) / PI) / 2.0 * n));

    // Clamp to valid range
    if (x < 0) x = 0;
    if (x >= static_cast<int>(n)) x = static_cast<int>(n) - 1;
    if (y < 0) y = 0;
    if (y >= static_cast<int>(n)) y = static_cast<int>(n) - 1;
}

// Convert longitude/latitude to quadbin cell at given resolution
inline uint64_t lonlat_to_cell(double lon, double lat, int z) {
    int x, y;
    lonlat_to_tile(lon, lat, z, x, y);
    return tile_to_cell(x, y, z);
}

// Convert tile coordinates to longitude/latitude (returns center of tile)
inline void tile_to_lonlat(int x, int y, int z, double &lon, double &lat) {
    double n = std::pow(2.0, z);
    lon = x / n * 360.0 - 180.0;
    double lat_rad = std::atan(std::sinh(PI * (1.0 - 2.0 * y / n)));
    lat = lat_rad * 180.0 / PI;
}

// Convert quadbin cell to longitude/latitude (center of cell)
inline void cell_to_lonlat(uint64_t cell, double &lon, double &lat) {
    int x, y, z;
    cell_to_tile(cell, x, y, z);
    // Return center of tile
    double n = std::pow(2.0, z);
    lon = (x + 0.5) / n * 360.0 - 180.0;
    double lat_rad = std::atan(std::sinh(PI * (1.0 - 2.0 * (y + 0.5) / n)));
    lat = lat_rad * 180.0 / PI;
}

// Get tile bounds in Web Mercator (EPSG:3857)
inline void tile_to_bbox_mercator(int x, int y, int z,
                                   double &min_x, double &min_y,
                                   double &max_x, double &max_y) {
    double n = std::pow(2.0, z);
    double tile_size = 2.0 * PI * EARTH_RADIUS / n;

    min_x = x * tile_size - PI * EARTH_RADIUS;
    max_x = (x + 1) * tile_size - PI * EARTH_RADIUS;

    // Y is flipped in Web Mercator tiles
    max_y = PI * EARTH_RADIUS - y * tile_size;
    min_y = PI * EARTH_RADIUS - (y + 1) * tile_size;
}

// Get tile bounds in WGS84 (EPSG:4326)
inline void tile_to_bbox_wgs84(int x, int y, int z,
                                double &min_lon, double &min_lat,
                                double &max_lon, double &max_lat) {
    double n = std::pow(2.0, z);

    min_lon = x / n * 360.0 - 180.0;
    max_lon = (x + 1) / n * 360.0 - 180.0;

    double min_lat_rad = std::atan(std::sinh(PI * (1.0 - 2.0 * (y + 1) / n)));
    double max_lat_rad = std::atan(std::sinh(PI * (1.0 - 2.0 * y / n)));

    min_lat = min_lat_rad * 180.0 / PI;
    max_lat = max_lat_rad * 180.0 / PI;
}

// Calculate pixel coordinates within a tile for a given lon/lat
inline void lonlat_to_pixel(double lon, double lat, int z, int tile_size,
                            int &pixel_x, int &pixel_y, int &tile_x, int &tile_y) {
    // Clamp latitude
    if (lat > MAX_LATITUDE) lat = MAX_LATITUDE;
    if (lat < -MAX_LATITUDE) lat = -MAX_LATITUDE;

    double n = std::pow(2.0, z);

    // Get fractional tile position
    double tile_x_frac = (lon + 180.0) / 360.0 * n;
    double lat_rad = lat * PI / 180.0;
    double tile_y_frac = (1.0 - std::asinh(std::tan(lat_rad)) / PI) / 2.0 * n;

    // Integer tile coordinates
    tile_x = static_cast<int>(std::floor(tile_x_frac));
    tile_y = static_cast<int>(std::floor(tile_y_frac));

    // Pixel within tile
    pixel_x = static_cast<int>((tile_x_frac - tile_x) * tile_size);
    pixel_y = static_cast<int>((tile_y_frac - tile_y) * tile_size);

    // Clamp to valid range
    if (pixel_x >= tile_size) pixel_x = tile_size - 1;
    if (pixel_y >= tile_size) pixel_y = tile_size - 1;
    if (pixel_x < 0) pixel_x = 0;
    if (pixel_y < 0) pixel_y = 0;
}

} // namespace quadbin
} // namespace duckdb
