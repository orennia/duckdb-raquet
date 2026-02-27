/// Quadbin geospatial index operations.
/// Ported from the C++ quadbin.hpp header.

const MAX_RESOLUTION: i32 = 26;
const EARTH_RADIUS: f64 = 6_378_137.0;
const MAX_LATITUDE: f64 = 85.051_128_779_806_6;
const PI: f64 = std::f64::consts::PI;

// QUADBIN encoding constants
const HEADER: u64 = 0x4000_0000_0000_0000;
const MODE: u64 = 0x0800_0000_0000_0000;
const FOOTER: u64 = 0x000F_FFFF_FFFF_FFFF;

// Magic numbers for Morton code interleaving
const B0: u64 = 0x5555_5555_5555_5555;
const B1: u64 = 0x3333_3333_3333_3333;
const B2: u64 = 0x0F0F_0F0F_0F0F_0F0F;
const B3: u64 = 0x00FF_00FF_00FF_00FF;
const B4: u64 = 0x0000_FFFF_0000_FFFF;
const B5: u64 = 0x0000_0000_FFFF_FFFF;

/// Extract resolution from a quadbin cell.
#[inline]
pub fn cell_to_resolution(cell: u64) -> i32 {
    ((cell >> 52) & 0x1F) as i32
}

/// Convert tile coordinates (x, y, z) to a quadbin cell.
pub fn tile_to_cell(x: i32, y: i32, z: i32) -> Result<u64, String> {
    if !(0..=MAX_RESOLUTION).contains(&z) {
        return Err(format!(
            "Resolution must be between 0 and {}",
            MAX_RESOLUTION
        ));
    }
    let mut ux = (x as u64) << (32 - z);
    let mut uy = (y as u64) << (32 - z);

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

    Ok(HEADER | MODE | ((z as u64) << 52) | ((ux | (uy << 1)) >> 12) | (FOOTER >> (z * 2)))
}

/// Convert a quadbin cell to tile coordinates (x, y, z).
pub fn cell_to_tile(cell: u64) -> (i32, i32, i32) {
    let z = cell_to_resolution(cell);
    let q = (cell & FOOTER) << 12;

    let mut ux = q;
    let mut uy = q >> 1;

    ux &= B0;
    uy &= B0;

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

    let x = (ux >> (32 - z)) as i32;
    let y = (uy >> (32 - z)) as i32;
    (x, y, z)
}

/// Convert lon/lat to tile coordinates at a given zoom level.
pub fn lonlat_to_tile(lon: f64, mut lat: f64, z: i32) -> (i32, i32) {
    if lat > MAX_LATITUDE {
        lat = MAX_LATITUDE;
    }
    if lat < -MAX_LATITUDE {
        lat = -MAX_LATITUDE;
    }
    let n = (2i64.pow(z as u32)) as f64;
    let mut x = ((lon + 180.0) / 360.0 * n).floor() as i32;
    let lat_rad = lat * PI / 180.0;
    let mut y = ((1.0 - lat_rad.tan().asinh() / PI) / 2.0 * n).floor() as i32;

    let max_coord = n as i32 - 1;
    if x < 0 {
        x = 0;
    }
    if x > max_coord {
        x = max_coord;
    }
    if y < 0 {
        y = 0;
    }
    if y > max_coord {
        y = max_coord;
    }
    (x, y)
}

/// Convert lon/lat to a quadbin cell at the given resolution.
pub fn lonlat_to_cell(lon: f64, lat: f64, z: i32) -> Result<u64, String> {
    let (x, y) = lonlat_to_tile(lon, lat, z);
    tile_to_cell(x, y, z)
}

/// Convert tile center to lon/lat.
pub fn tile_to_lonlat(x: i32, y: i32, z: i32) -> (f64, f64) {
    let n = (2i64.pow(z as u32)) as f64;
    let lon = x as f64 / n * 360.0 - 180.0;
    let lat_rad = (PI * (1.0 - 2.0 * y as f64 / n)).sinh().atan();
    (lon, lat_rad * 180.0 / PI)
}

/// Convert a quadbin cell to lon/lat (center of cell).
pub fn cell_to_lonlat(cell: u64) -> (f64, f64) {
    let (x, y, z) = cell_to_tile(cell);
    let n = (2i64.pow(z as u32)) as f64;
    let lon = (x as f64 + 0.5) / n * 360.0 - 180.0;
    let lat_rad = (PI * (1.0 - 2.0 * (y as f64 + 0.5) / n)).sinh().atan();
    (lon, lat_rad * 180.0 / PI)
}

/// Get tile bounds in WGS84 (EPSG:4326).
pub fn tile_to_bbox_wgs84(x: i32, y: i32, z: i32) -> (f64, f64, f64, f64) {
    let n = (2i64.pow(z as u32)) as f64;
    let min_lon = x as f64 / n * 360.0 - 180.0;
    let max_lon = (x as f64 + 1.0) / n * 360.0 - 180.0;
    let min_lat_rad = (PI * (1.0 - 2.0 * (y as f64 + 1.0) / n)).sinh().atan();
    let max_lat_rad = (PI * (1.0 - 2.0 * y as f64 / n)).sinh().atan();
    (
        min_lon,
        min_lat_rad * 180.0 / PI,
        max_lon,
        max_lat_rad * 180.0 / PI,
    )
}

/// Get tile bounds in Web Mercator (EPSG:3857).
pub fn tile_to_bbox_mercator(x: i32, y: i32, z: i32) -> (f64, f64, f64, f64) {
    let n = (2i64.pow(z as u32)) as f64;
    let tile_size = 2.0 * PI * EARTH_RADIUS / n;
    let min_x = x as f64 * tile_size - PI * EARTH_RADIUS;
    let max_x = (x as f64 + 1.0) * tile_size - PI * EARTH_RADIUS;
    let max_y = PI * EARTH_RADIUS - y as f64 * tile_size;
    let min_y = PI * EARTH_RADIUS - (y as f64 + 1.0) * tile_size;
    (min_x, min_y, max_x, max_y)
}

/// Get parent cell at the given resolution.
pub fn cell_to_parent_at(cell: u64, parent_resolution: i32) -> Result<u64, String> {
    let current_res = cell_to_resolution(cell);
    if parent_resolution < 0 || parent_resolution > current_res {
        return Err(format!(
            "Parent resolution must be between 0 and {}",
            current_res
        ));
    }
    if parent_resolution == current_res {
        return Ok(cell);
    }
    let (x, y, z) = cell_to_tile(cell);
    let shift = z - parent_resolution;
    tile_to_cell(x >> shift, y >> shift, parent_resolution)
}

/// Get parent cell at resolution - 1.
pub fn cell_to_parent(cell: u64) -> Result<u64, String> {
    let current_res = cell_to_resolution(cell);
    if current_res == 0 {
        return Ok(cell);
    }
    cell_to_parent_at(cell, current_res - 1)
}

/// Get all child cells at the given resolution.
pub fn cell_to_children_at(cell: u64, child_resolution: i32) -> Result<Vec<u64>, String> {
    let current_res = cell_to_resolution(cell);
    if child_resolution <= current_res || child_resolution > MAX_RESOLUTION {
        return Err(format!(
            "Child resolution must be greater than current ({}) and <= {}",
            current_res, MAX_RESOLUTION
        ));
    }
    let (x, y, _z) = cell_to_tile(cell);
    let res_diff = child_resolution - current_res;
    let children_per_dim = 1i32 << res_diff;
    let base_x = x << res_diff;
    let base_y = y << res_diff;
    let mut children = Vec::with_capacity((children_per_dim * children_per_dim) as usize);
    for dy in 0..children_per_dim {
        for dx in 0..children_per_dim {
            children.push(tile_to_cell(base_x + dx, base_y + dy, child_resolution)?);
        }
    }
    Ok(children)
}

/// Get the 4 immediate children of a cell (at resolution + 1).
pub fn cell_to_children(cell: u64) -> Result<Vec<u64>, String> {
    let res = cell_to_resolution(cell);
    cell_to_children_at(cell, res + 1)
}

/// Get k-ring neighbors (cells within grid distance k).
pub fn cell_kring(cell: u64, k: i32) -> Result<Vec<u64>, String> {
    if k < 0 {
        return Err("K must be non-negative".to_string());
    }
    let (x, y, z) = cell_to_tile(cell);
    let max_coord = (1i32 << z) - 1;
    let mut neighbors = Vec::new();
    for dy in -k..=k {
        for dx in -k..=k {
            let nx = x + dx;
            let ny = y + dy;
            if nx < 0 || nx > max_coord || ny < 0 || ny > max_coord {
                continue;
            }
            neighbors.push(tile_to_cell(nx, ny, z)?);
        }
    }
    Ok(neighbors)
}

/// Get the 4 sibling cells (other children of the same parent).
pub fn cell_siblings(cell: u64) -> Result<Vec<u64>, String> {
    let z = cell_to_resolution(cell);
    if z == 0 {
        return Ok(vec![cell; 4]);
    }
    let parent = cell_to_parent(cell)?;
    cell_to_children(parent)
}

/// Calculate pixel coordinates within a tile for a given lon/lat.
pub fn lonlat_to_pixel(lon: f64, mut lat: f64, z: i32, tile_size: i32) -> (i32, i32, i32, i32) {
    if lat > MAX_LATITUDE {
        lat = MAX_LATITUDE;
    }
    if lat < -MAX_LATITUDE {
        lat = -MAX_LATITUDE;
    }
    let n = (2i64.pow(z as u32)) as f64;
    let tile_x_frac = (lon + 180.0) / 360.0 * n;
    let lat_rad = lat * PI / 180.0;
    let tile_y_frac = (1.0 - lat_rad.tan().asinh() / PI) / 2.0 * n;

    let tile_x = tile_x_frac.floor() as i32;
    let tile_y = tile_y_frac.floor() as i32;

    let mut pixel_x = ((tile_x_frac - tile_x as f64) * tile_size as f64) as i32;
    let mut pixel_y = ((tile_y_frac - tile_y as f64) * tile_size as f64) as i32;

    if pixel_x >= tile_size {
        pixel_x = tile_size - 1;
    }
    if pixel_y >= tile_size {
        pixel_y = tile_size - 1;
    }
    if pixel_x < 0 {
        pixel_x = 0;
    }
    if pixel_y < 0 {
        pixel_y = 0;
    }

    (pixel_x, pixel_y, tile_x, tile_y)
}

// ============================================================================
// Polyfill: generate all quadbin cells covering a WKB geometry bounding box
// ============================================================================

/// Check if a point is inside a polygon ring (ray casting).
pub fn point_in_ring(px: f64, py: f64, coords: &[f64]) -> bool {
    let n = coords.len() / 2;
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let xi = coords[i * 2];
        let yi = coords[i * 2 + 1];
        let xj = coords[j * 2];
        let yj = coords[j * 2 + 1];
        if ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Parse a WKB polygon ring and return (coords, bytes_consumed).
fn parse_wkb_ring(data: &[u8], little_endian: bool) -> Option<(Vec<f64>, usize)> {
    if data.len() < 4 {
        return None;
    }
    let num_points = read_u32(data, little_endian)? as usize;
    let needed = 4 + num_points * 16;
    if data.len() < needed {
        return None;
    }
    let mut coords = Vec::with_capacity(num_points * 2);
    for i in 0..num_points {
        let off = 4 + i * 16;
        let x = read_f64(&data[off..], little_endian)?;
        let y = read_f64(&data[off + 8..], little_endian)?;
        coords.push(x);
        coords.push(y);
    }
    Some((coords, needed))
}

pub fn read_u32(data: &[u8], little_endian: bool) -> Option<u32> {
    if data.len() < 4 {
        return None;
    }
    let b = [data[0], data[1], data[2], data[3]];
    Some(if little_endian {
        u32::from_le_bytes(b)
    } else {
        u32::from_be_bytes(b)
    })
}

pub fn read_f64(data: &[u8], little_endian: bool) -> Option<f64> {
    if data.len() < 8 {
        return None;
    }
    let b: [u8; 8] = data[..8].try_into().ok()?;
    Some(if little_endian {
        f64::from_le_bytes(b)
    } else {
        f64::from_be_bytes(b)
    })
}

/// Check if a point is inside a WKB polygon (with holes).
fn point_in_wkb_polygon(px: f64, py: f64, data: &[u8], little_endian: bool) -> Option<bool> {
    if data.len() < 4 {
        return Some(false);
    }
    let num_rings = read_u32(data, little_endian)? as usize;
    if num_rings == 0 {
        return Some(false);
    }
    let mut offset = 4usize;
    let (outer, consumed) = parse_wkb_ring(&data[offset..], little_endian)?;
    offset += consumed;
    if !point_in_ring(px, py, &outer) {
        return Some(false);
    }
    for _ in 1..num_rings {
        let (hole, consumed) = parse_wkb_ring(&data[offset..], little_endian)?;
        offset += consumed;
        if point_in_ring(px, py, &hole) {
            return Some(false);
        }
    }
    Some(true)
}

/// Parse bounding box and geometry type from a WKB geometry.
/// Returns (min_x, min_y, max_x, max_y, geom_data_for_point_test, little_endian, is_polygon).
#[allow(clippy::type_complexity)]
pub fn parse_wkb_bbox(data: &[u8]) -> Option<(f64, f64, f64, f64, Vec<u8>, bool, bool)> {
    if data.len() < 5 {
        return None;
    }
    let byte_order = data[0];
    if byte_order > 1 {
        return None;
    }
    let little_endian = byte_order == 1;
    let geom_type = read_u32(&data[1..], little_endian)?;
    let base_type = geom_type & 0xFF;
    let mut offset = 5usize;
    // Skip SRID if present
    if geom_type & 0x2000_0000 != 0 {
        offset += 4;
    }

    let mut min_x = f64::MAX;
    let mut min_y = f64::MAX;
    let mut max_x = f64::MIN;
    let mut max_y = f64::MIN;

    match base_type {
        3 => {
            // POLYGON: read outer ring for bbox
            if offset + 4 > data.len() {
                return None;
            }
            let num_rings = read_u32(&data[offset..], little_endian)? as usize;
            offset += 4;
            for r in 0..num_rings {
                if offset + 4 > data.len() {
                    return None;
                }
                let num_pts = read_u32(&data[offset..], little_endian)? as usize;
                offset += 4;
                for _ in 0..num_pts {
                    if offset + 16 > data.len() {
                        return None;
                    }
                    let x = read_f64(&data[offset..], little_endian)?;
                    let y = read_f64(&data[offset + 8..], little_endian)?;
                    if r == 0 {
                        if x < min_x {
                            min_x = x;
                        }
                        if x > max_x {
                            max_x = x;
                        }
                        if y < min_y {
                            min_y = y;
                        }
                        if y > max_y {
                            max_y = y;
                        }
                    }
                    offset += 16;
                }
            }
        }
        6 => {
            // MULTIPOLYGON
            if offset + 4 > data.len() {
                return None;
            }
            let num_polys = read_u32(&data[offset..], little_endian)? as usize;
            offset += 4;
            for _ in 0..num_polys {
                if offset + 5 > data.len() {
                    return None;
                }
                let poly_le = data[offset] == 1;
                offset += 5; // byte_order + geom_type
                if offset + 4 > data.len() {
                    return None;
                }
                let num_rings = read_u32(&data[offset..], poly_le)? as usize;
                offset += 4;
                for r in 0..num_rings {
                    if offset + 4 > data.len() {
                        return None;
                    }
                    let num_pts = read_u32(&data[offset..], poly_le)? as usize;
                    offset += 4;
                    for _ in 0..num_pts {
                        if offset + 16 > data.len() {
                            return None;
                        }
                        let x = read_f64(&data[offset..], poly_le)?;
                        let y = read_f64(&data[offset + 8..], poly_le)?;
                        if r == 0 {
                            if x < min_x {
                                min_x = x;
                            }
                            if x > max_x {
                                max_x = x;
                            }
                            if y < min_y {
                                min_y = y;
                            }
                            if y > max_y {
                                max_y = y;
                            }
                        }
                        offset += 16;
                    }
                }
            }
        }
        _ => return None,
    }

    if min_x > max_x {
        return None;
    }
    let geom_data = data.to_vec();
    Some((min_x, min_y, max_x, max_y, geom_data, little_endian, true))
}

/// Polyfill: return all quadbin cells at `resolution` that intersect the WKB geometry.
pub fn quadbin_polyfill(wkb: &[u8], resolution: i32) -> Result<Vec<u64>, String> {
    if !(0..=MAX_RESOLUTION).contains(&resolution) {
        return Err(format!(
            "Resolution must be between 0 and {}",
            MAX_RESOLUTION
        ));
    }
    let parsed = match parse_wkb_bbox(wkb) {
        Some(p) => p,
        None => return Ok(vec![]),
    };
    let (min_lon, min_lat, max_lon, max_lat, geom_data, little_endian, _) = parsed;

    // Get tile range covering the bbox
    let (x_min, y_max) = lonlat_to_tile(min_lon, min_lat, resolution);
    let (x_max, y_min) = lonlat_to_tile(max_lon, max_lat, resolution);

    let mut result = Vec::new();
    for ty in y_min..=y_max {
        for tx in x_min..=x_max {
            // Get center of this tile
            let (cx, cy) = tile_center_lonlat(tx, ty, resolution);
            // Check if center is inside polygon
            if point_in_wkb_geometry(cx, cy, &geom_data, little_endian) {
                if let Ok(cell) = tile_to_cell(tx, ty, resolution) {
                    result.push(cell);
                }
            }
        }
    }
    Ok(result)
}

fn tile_center_lonlat(x: i32, y: i32, z: i32) -> (f64, f64) {
    let n = (2i64.pow(z as u32)) as f64;
    let lon = (x as f64 + 0.5) / n * 360.0 - 180.0;
    let lat_rad = (PI * (1.0 - 2.0 * (y as f64 + 0.5) / n)).sinh().atan();
    (lon, lat_rad * 180.0 / PI)
}

fn point_in_wkb_geometry(px: f64, py: f64, data: &[u8], _default_le: bool) -> bool {
    if data.len() < 5 {
        return false;
    }
    let le = data[0] == 1;
    let geom_type = match read_u32(&data[1..], le) {
        Some(t) => t,
        None => return false,
    };
    let base_type = geom_type & 0xFF;
    let mut offset = 5usize;
    if geom_type & 0x2000_0000 != 0 {
        offset += 4;
    }

    match base_type {
        3 => {
            // POLYGON
            point_in_wkb_polygon(px, py, &data[offset..], le).unwrap_or(false)
        }
        6 => {
            // MULTIPOLYGON
            if offset + 4 > data.len() {
                return false;
            }
            let num_polys = match read_u32(&data[offset..], le) {
                Some(n) => n as usize,
                None => return false,
            };
            offset += 4;
            for _ in 0..num_polys {
                if offset + 5 > data.len() {
                    break;
                }
                let poly_le = data[offset] == 1;
                offset += 5;
                if let Some(inside) = point_in_wkb_polygon(px, py, &data[offset..], poly_le) {
                    if inside {
                        return true;
                    }
                    // Advance past this polygon
                    if offset + 4 > data.len() {
                        break;
                    }
                    let num_rings = match read_u32(&data[offset..], poly_le) {
                        Some(n) => n as usize,
                        None => break,
                    };
                    offset += 4;
                    for _ in 0..num_rings {
                        if offset + 4 > data.len() {
                            break;
                        }
                        let num_pts = match read_u32(&data[offset..], poly_le) {
                            Some(n) => n as usize,
                            None => break,
                        };
                        offset += 4 + num_pts * 16;
                    }
                }
            }
            false
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tile_to_cell_root() {
        let cell = tile_to_cell(0, 0, 0).unwrap();
        assert_eq!(cell, 5_192_650_370_358_181_887);
    }

    #[test]
    fn test_roundtrip() {
        for z in 0..=13 {
            let x = 100 % (1 << z);
            let y = 50 % (1 << z);
            let cell = tile_to_cell(x, y, z).unwrap();
            let (rx, ry, rz) = cell_to_tile(cell);
            assert_eq!((x, y, z), (rx, ry, rz));
        }
    }

    #[test]
    fn test_resolution() {
        let cell = tile_to_cell(0, 0, 5).unwrap();
        assert_eq!(cell_to_resolution(cell), 5);
    }
}
