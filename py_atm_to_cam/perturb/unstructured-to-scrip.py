import numpy as np
from netCDF4 import Dataset
import xarray as xr
from scipy.spatial import SphericalVoronoi
import matplotlib.pyplot as plt
import warnings

def lat_lon_to_xyz(lat, lon):
    """Convert latitude/longitude to 3D Cartesian coordinates on unit sphere."""
    lat_rad, lon_rad = np.radians(lat), np.radians(lon)
    cos_lat = np.cos(lat_rad)
    return np.column_stack((
        cos_lat * np.cos(lon_rad),
        cos_lat * np.sin(lon_rad),
        np.sin(lat_rad)
    ))

def xyz_to_lat_lon(xyz):
    """Convert 3D Cartesian coordinates to latitude/longitude."""
    xyz_normalized = xyz / np.linalg.norm(xyz, axis=1)[:, np.newaxis]
    lat = np.degrees(np.arcsin(xyz_normalized[:, 2]))
    lon = np.degrees(np.arctan2(xyz_normalized[:, 1], xyz_normalized[:, 0]))
    return lat, lon

def sort_vertices_spherical(vertices, center_point):
    """
    Sort vertices counterclockwise around a center point on a sphere.

    Parameters:
    -----------
    vertices : ndarray
        Array of shape (N, 3) containing vertex coordinates in 3D
    center_point : ndarray
        Array of shape (3,) containing the center point coordinates

    Returns:
    --------
    ndarray
        Indices that will sort vertices counterclockwise
    """
    # Normalize all points to ensure they're on unit sphere
    center_norm = center_point / np.linalg.norm(center_point)
    vertices_norm = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis]

    # Get North pole vector (or different reference if too close to poles)
    if abs(np.dot(center_norm, [0, 0, 1])) > 0.99:
        ref = np.array([1, 0, 0])
    else:
        ref = np.array([0, 0, 1])

    # Get reference direction (tangent to sphere at center_point)
    ref_dir = ref - center_norm * np.dot(ref, center_norm)
    ref_dir = ref_dir / np.linalg.norm(ref_dir)

    # Get vector perpendicular to both center_norm and ref_dir
    perp_dir = np.cross(center_norm, ref_dir)

    # For each vertex, compute the tangent vector from center to vertex
    # (removing the radial component)
    tangent_vectors = vertices_norm - center_norm * \
                     np.dot(vertices_norm, center_norm)[:, np.newaxis]

    # Normalize tangent vectors
    norms = np.linalg.norm(tangent_vectors, axis=1)
    valid = norms > 1e-10
    tangent_vectors[valid] /= norms[valid, np.newaxis]

    # Compute angles using atan2
    # Project onto our reference directions
    x_proj = np.dot(tangent_vectors, ref_dir)
    y_proj = np.dot(tangent_vectors, perp_dir)
    angles = np.arctan2(y_proj, x_proj)

    # Handle numerical instability
    angles[~valid] = 0

    return np.argsort(angles)

def normalize_lon_0_360(lon):
    """
    Normalize longitude to [0,360) range.
    """
    return lon % 360

def find_spherical_corners(lat, lon, debug=False):
    """
    Find corners for each point using spherical Voronoi tessellation with
    guaranteed counterclockwise ordering.
    """
    # Normalize input longitudes first
    lon = normalize_lon_0_360(lon)

    # Convert to 3D coordinates
    points_3d = lat_lon_to_xyz(lat, lon)

    # Create spherical Voronoi tessellation
    sv = SphericalVoronoi(points_3d)

    grid_size = len(lat)

    # Determine maximum number of corners across all regions
    max_corners = max(len(region) for region in sv.regions)
    if debug:
        print(f"Maximum number of corners in any cell: {max_corners}")

    corner_lats = np.full((grid_size, max_corners), np.nan)
    corner_lons = np.full((grid_size, max_corners), np.nan)

    # Get regions and vertices
    regions = sv.regions
    vertices = sv.vertices

    if debug:
        print(f"Found {len(regions)} regions")
        print(f"Found {len(vertices)} vertices")

    # Process each point and its corresponding region
    for i in range(grid_size):
        if debug and i % 1000 == 0:
            print(f"Processing point {i}/{grid_size}")

        region = regions[i]
        if len(region) > 0:
            # Get vertices for this region
            region_vertices = vertices[region]

            # Sort vertices counterclockwise
            sorted_idx = sort_vertices_spherical(region_vertices, points_3d[i])
            sorted_vertices = region_vertices[sorted_idx]

            # Convert to lat/lon
            rlat, rlon = xyz_to_lat_lon(sorted_vertices)

            # Normalize longitudes to [0,360)
            rlon = normalize_lon_0_360(rlon)

            # Handle cases where corners cross the prime meridian
            if lon[i] < 30 or lon[i] > 330:
                mean_lon = np.mean(rlon)
                if abs(mean_lon - lon[i]) > 180:
                    if lon[i] < 30:
                        rlon[rlon > 180] -= 360
                    else:
                        rlon[rlon < 180] += 360
                rlon = normalize_lon_0_360(rlon)

            # Store all corners
            n_corners = len(rlat)
            corner_lats[i, :n_corners] = rlat
            corner_lons[i, :n_corners] = rlon

            # Fill remaining slots with the last corner
            if n_corners < max_corners:
                corner_lats[i, n_corners:] = rlat[-1]
                corner_lons[i, n_corners:] = rlon[-1]

    return corner_lats, corner_lons

def calculate_spherical_areas(corner_lats, corner_lons):
    """Calculate areas of spherical polygons more efficiently."""
    grid_size = len(corner_lats)
    areas = np.zeros(grid_size)

    # Convert all corners to xyz at once
    xyz = lat_lon_to_xyz(corner_lats.ravel(), corner_lons.ravel())
    xyz = xyz.reshape(grid_size, -1, 3)  # Reshape back to (grid_size, corners, 3)

    # Calculate areas using vectorized operations where possible
    for i in range(grid_size):
        valid = ~np.isnan(corner_lats[i])
        if np.sum(valid) >= 3:
            points = xyz[i, valid]
            area = 0
            for j in range(len(points)-1):
                area += np.arctan2(
                    np.abs(np.dot(points[0], np.cross(points[j], points[j+1]))),
                    1 + np.dot(points[0], points[j]) + np.dot(points[0], points[j+1]) +
                    np.dot(points[j], points[j+1])
                )
            areas[i] = abs(area)

    return areas

def unstructured_to_scrip(filename, lat, lon, grid_corners=None, grid_mask=None,
                         grid_area=None, title=None, debug=False):
    """
    Convert latitude and longitude arrays to SCRIP unstructured grid format.
    """
    # Convert inputs to numpy arrays and check dimensions
    lat = np.asarray(lat).ravel()
    lon = np.asarray(lon).ravel()

    if lat.size != lon.size:
        raise ValueError("Latitude and longitude must have the same number of elements")

    grid_size = lat.size
    grid_rank = 1  # Unstructured grid

    # Create corners if not provided
    if grid_corners is None:
        if debug:
            print("Creating spherical Voronoi tessellation...")
        corner_lats, corner_lons = find_spherical_corners(lat, lon, debug)
        grid_corners = (corner_lats, corner_lons)

    corner_lats, corner_lons = grid_corners
    num_corners = corner_lats.shape[1]

    # Calculate areas if not provided
    if grid_area is None:
        if debug:
            print("Calculating spherical areas...")
        grid_area = calculate_spherical_areas(corner_lats, corner_lons)

    # Create default mask if not provided
    if grid_mask is None:
        grid_mask = np.ones(grid_size, dtype=np.int32)

    if debug:
        print("Writing SCRIP file...")

    # Create the NetCDF file
    with Dataset(filename, 'w', format='NETCDF4') as nc:
        # Create dimensions
        nc.createDimension('grid_size', grid_size)
        nc.createDimension('grid_corners', num_corners)
        nc.createDimension('grid_rank', grid_rank)

        # Create variables
        nc.createVariable('grid_dims', 'i4', ('grid_rank',))

        # Center coordinates
        center_lat = nc.createVariable('grid_center_lat', 'f8', ('grid_size',))
        center_lat.units = 'degrees'
        center_lat[:] = lat

        center_lon = nc.createVariable('grid_center_lon', 'f8', ('grid_size',))
        center_lon.units = 'degrees'
        center_lon[:] = lon

        # Corner coordinates
        corner_lat = nc.createVariable('grid_corner_lat', 'f8',
                                     ('grid_size', 'grid_corners'),
                                     fill_value=-9999.)
        corner_lat.units = 'degrees'
        corner_lat[:] = corner_lats

        corner_lon = nc.createVariable('grid_corner_lon', 'f8',
                                     ('grid_size', 'grid_corners'),
                                     fill_value=-9999.)
        corner_lon.units = 'degrees'
        corner_lon[:] = corner_lons

        # Mask
        mask = nc.createVariable('grid_imask', 'i4', ('grid_size',),
                               fill_value=-9999)
        mask.units = 'unitless'
        mask[:] = grid_mask

        # Grid dimensions
        nc.variables['grid_dims'][:] = [grid_size]

        # Area
        area = nc.createVariable('grid_area', 'f8', ('grid_size',))
        area.units = 'radians^2'
        area.long_name = 'area weights'
        area[:] = grid_area

        # Global attributes
        nc.Conventions = 'SCRIP'
        nc.title = title if title else f'SCRIP Grid ({grid_size} points)'
        nc.created_by = 'Python unstructured_to_scrip function'

        if debug:
            print(f"Created SCRIP file with {grid_size} points and {num_corners} corners per cell")

def spherical_polygon_area(lats, lons):
    """
    Calculate the area of a spherical polygon using L'Huilier's formula.

    Parameters:
    -----------
    lats, lons : array-like
        Arrays of latitude and longitude in radians

    Returns:
    -------
    float
        Area in steradians
    """
    if len(lats) < 3:
        return 0.0

    # Convert to 3D coordinates
    points = np.column_stack([
        np.cos(lats) * np.cos(lons),
        np.cos(lats) * np.sin(lons),
        np.sin(lats)
    ])

    # Calculate the area using spherical excess
    total = 0
    n = len(points)
    for i in range(n):
        j = (i + 1) % n
        k = (i + 2) % n

        # Get vectors for triangle
        a = points[i]
        b = points[j]
        c = points[k]

        # Calculate the spherical excess for this triangle
        numerator = abs(np.dot(a, np.cross(b, c)))
        denominator = (1 + np.dot(a, b) + np.dot(b, c) + np.dot(c, a))
        angle = 2 * np.arctan2(numerator, denominator)
        total += angle

    # Subtract out extra triangles
    area = total - (n - 2) * np.pi
    return abs(area)

def check_sphere_coverage(lat, lon, corner_lats, corner_lons, debug=False, create_plot=False):
    """
    Check if the cells properly cover the sphere.
    """
    # Convert to radians
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    corner_lats_rad = np.radians(corner_lats)
    corner_lons_rad = np.radians(corner_lons)

    # Calculate areas
    areas = np.zeros(len(lat))
    for i in range(len(lat)):
        valid = ~np.isnan(corner_lats_rad[i])
        if np.sum(valid) >= 3:
            areas[i] = spherical_polygon_area(
                corner_lats_rad[i, valid],
                corner_lons_rad[i, valid]
            )

    # Rest of the function remains the same...
    total_area = np.sum(areas)
    sphere_area = 4 * np.pi
    coverage_fraction = total_area / sphere_area

    print(f"DEBUG: Total computed area: {total_area:.4f}")
    print(f"DEBUG: Sphere area: {sphere_area:.4f}")
    print(f"DEBUG: Coverage fraction: {coverage_fraction:.4f}")
    print(f"DEBUG: Number of non-zero areas: {np.sum(areas > 0)}")
    print(f"DEBUG: Area range: [{np.min(areas):.6f}, {np.max(areas):.6f}]")

    # Check for gaps between cells
    has_gaps = coverage_fraction < 0.99  # Allow 1% tolerance

    # Check for cell validity
    valid_cells = np.all(~np.isnan(corner_lats), axis=1) & np.all(~np.isnan(corner_lons), axis=1)
    invalid_cells = np.where(~valid_cells)[0]

    # Check for suspicious areas
    mean_area = np.mean(areas[areas > 0])
    std_area = np.std(areas[areas > 0])
    suspicious_areas = np.where((areas > mean_area + 3*std_area) |
                              (areas < mean_area - 3*std_area))[0]

    if debug and len(suspicious_areas) > 0:
        print("\nSuspicious Cells Found:")
        print("----------------------")
        print(f"Mean area: {mean_area:.6f}")
        print(f"Std dev:  {std_area:.6f}")
        print(f"Expected range: [{mean_area - 3*std_area:.6f}, {mean_area + 3*std_area:.6f}]")
        print("\nDetailed cell information:")
        print("  Cell ID    Latitude    Longitude        Area      Deviation")
        print("--------------------------------------------------------")
        for cell_id in suspicious_areas:
            deviation = (areas[cell_id] - mean_area) / std_area
            print(f"  {cell_id:6d}    {lat[cell_id]:9.4f}    {lon[cell_id]:9.4f}    {areas[cell_id]:9.6f}    {deviation:9.2f}σ")

    if create_plot:
        # Create visualization
        plt.figure(figsize=(15, 5))

        # Plot cell areas
        plt.subplot(131)
        plt.hist(areas, bins=50)
        plt.title('Cell Areas Distribution')
        plt.xlabel('Area')
        plt.ylabel('Count')

        # Plot points on sphere
        plt.subplot(132)
        plt.scatter(lon, lat, c=areas, cmap='viridis')
        plt.colorbar(label='Area')
        plt.title('Cell Areas on Sphere')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        # Plot coverage gaps if any
        plt.subplot(133)
        plt.scatter(lon, lat, c='blue', alpha=0.5, label='Valid cells')
        if len(invalid_cells) > 0:
            plt.scatter(lon[invalid_cells], lat[invalid_cells],
                       c='red', label='Invalid cells')
        if len(suspicious_areas) > 0:
            plt.scatter(lon[suspicious_areas], lat[suspicious_areas],
                       c='yellow', label='Suspicious areas')
        plt.legend()
        plt.title('Coverage Issues')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        plt.tight_layout()
        plt.show()

    return {
        'total_area': total_area,
        'coverage_fraction': coverage_fraction,
        'has_gaps': has_gaps,
        'num_invalid_cells': len(invalid_cells),
        'invalid_cell_indices': invalid_cells,
        'suspicious_areas': suspicious_areas,
        'mean_cell_area': mean_area,
        'std_cell_area': std_area,
        'min_area': np.min(areas),
        'max_area': np.max(areas)
    }

def test_grid_coverage(filename):
    """
    Test the coverage of a SCRIP grid file.

    Parameters:
    -----------
    filename : str
        Path to SCRIP grid file
    """
    with Dataset(filename, 'r') as nc:
        lat = nc.variables['grid_center_lat'][:]
        lon = nc.variables['grid_center_lon'][:]
        corner_lats = nc.variables['grid_corner_lat'][:]
        corner_lons = nc.variables['grid_corner_lon'][:]

    results = check_sphere_coverage(lat, lon, corner_lats, corner_lons, debug=True)

    print("\nGrid Coverage Analysis:")
    print(f"Total area covered: {results['total_area']:.2f} (should be close to 4π ≈ 12.57)")
    print(f"Coverage fraction: {results['coverage_fraction']:.2%}")
    print(f"Has gaps: {results['has_gaps']}")
    print(f"Number of invalid cells: {results['num_invalid_cells']}")
    if results['num_invalid_cells'] > 0:
        print(f"Invalid cell indices: {results['invalid_cell_indices']}")
    print(f"Number of suspicious areas: {len(results['suspicious_areas'])}")
    print(f"\nCell area statistics:")
    print(f"Mean area: {results['mean_cell_area']:.6f}")
    print(f"Std dev: {results['std_cell_area']:.6f}")
    print(f"Min area: {results['min_area']:.6f}")
    print(f"Max area: {results['max_area']:.6f}")

    return results

def debug_cell(cell_id, corner_lats, corner_lons, center_lat, center_lon):
    """
    Debug a specific cell's geometry and ordering.

    Parameters:
    -----------
    cell_id : int
        The ID of the cell to debug
    corner_lats, corner_lons : ndarray
        Arrays containing corner coordinates
    center_lat, center_lon : float
        Center coordinates of the cell
    """
    import matplotlib.pyplot as plt

    # Get the corners for this cell
    clat = corner_lats[cell_id]
    clon = corner_lons[cell_id]

    # Remove any NaN values
    valid = ~np.isnan(clat)
    clat = clat[valid]
    clon = clon[valid]

    print(f"\nDEBUG INFO FOR CELL {cell_id}")
    print("----------------------------------------")
    print(f"Center coordinates: lat={center_lat:.6f}, lon={center_lon:.6f}")
    print("\nCorner coordinates:")
    for i, (lat, lon) in enumerate(zip(clat, clon)):
        print(f"Corner {i}: lat={lat:.6f}, lon={lon:.6f}")

    # Calculate angles between successive corners
    center_xyz = lat_lon_to_xyz(np.array([center_lat]), np.array([center_lon]))[0]
    corner_xyz = lat_lon_to_xyz(clat, clon)

    angles = []
    for i in range(len(corner_xyz)):
        v1 = corner_xyz[i]
        v2 = corner_xyz[(i+1) % len(corner_xyz)]
        # Compute signed angle between vectors in the tangent plane
        angle = np.arctan2(
            np.dot(np.cross(v1-center_xyz, v2-center_xyz), center_xyz),
            np.dot(v1-center_xyz, v2-center_xyz)
        )
        angles.append(np.degrees(angle))

    print("\nAngles between successive corners (degrees):")
    angle_sum = 0
    for i, angle in enumerate(angles):
        angle_sum += angle
        print(f"Angle {i}->{i+1}: {angle:.2f}")
    print(f"Sum of angles: {angle_sum:.2f} (should be close to 360 for CCW)")

#     # Plot the cell
#     plt.figure(figsize=(10, 10))
#
#     # Plot corners
#     plt.scatter(clon, clat, c='blue', label='Corners')
#
#     # Plot center
#     plt.scatter(center_lon, center_lat, c='red', label='Center')
#
#     # Connect corners with lines
#     for i in range(len(clat)):
#         next_i = (i + 1) % len(clat)
#         plt.plot([clon[i], clon[next_i]], [clat[i], clat[next_i]], 'b-')
#         # Add corner numbers
#         plt.annotate(f'C{i}', (clon[i], clat[i]))
#
#     plt.title(f'Cell {cell_id} Geometry')
#     plt.xlabel('Longitude')
#     plt.ylabel('Latitude')
#     plt.legend()
#     plt.grid(True)
#     plt.show()

    # Check for potential issues
    print("\nPotential issues:")

    # Check for very close corners
    for i in range(len(corner_xyz)):
        for j in range(i+1, len(corner_xyz)):
            dist = np.linalg.norm(corner_xyz[i] - corner_xyz[j])
            if dist < 1e-6:
                print(f"WARNING: Corners {i} and {j} are very close together")

    # Check for colinear points
    for i in range(len(corner_xyz)):
        v1 = corner_xyz[i] - center_xyz
        v2 = corner_xyz[(i+1) % len(corner_xyz)] - center_xyz
        cross_prod = np.linalg.norm(np.cross(v1, v2))
        if cross_prod < 1e-6:
            print(f"WARNING: Corners {i} and {i+1} are nearly colinear with center")

    # Check for clockwise ordering
    if angle_sum < 0:
        print("WARNING: Corners appear to be ordered clockwise")

    return angles

def test_specific_cell(filename, cell_id):
    """
    Test a specific cell in a SCRIP grid file.
    """
    with Dataset(filename, 'r') as nc:
        corner_lats = nc.variables['grid_corner_lat'][:]
        corner_lons = nc.variables['grid_corner_lon'][:]
        center_lat = nc.variables['grid_center_lat'][cell_id]
        center_lon = nc.variables['grid_center_lon'][cell_id]

    debug_cell(cell_id, corner_lats, corner_lons, center_lat, center_lon)

def validate_scrip_file(filename):
    """
    Validate a SCRIP grid file and print results.

    Parameters:
    -----------
    filename : str
        Path to SCRIP file to validate
    """
    print("\nValidating SCRIP file:", filename)
    print("============================")

    # Get coverage results
    results = test_grid_coverage(filename)

    # Check critical issues
    has_issues = False

    if results['coverage_fraction'] < 0.99:
        print(f"ERROR: Incomplete coverage: {results['coverage_fraction']:.1%}")
        has_issues = True

    if results['num_invalid_cells'] > 0:
        print(f"ERROR: Found {results['num_invalid_cells']} invalid cells")
        print(f"       First few invalid cells: {results['invalid_cell_indices'][:5]}")
        has_issues = True

    if len(results['suspicious_areas']) > 0:
        print(f"WARNING: Found {len(results['suspicious_areas'])} cells with suspicious areas")

    # Print basic statistics
    print(f"\nGrid size: {results['total_area']:.2f}")
    print(f"Mean cell area: {results['mean_cell_area']:.6f}")
    print(f"Area range: [{results['min_area']:.6f}, {results['max_area']:.6f}]")

    if not has_issues:
        print("\nNo critical issues found.")

    return not has_issues  # Return True if valid, False if issues found

if __name__ == '__main__':
    # Example usage

    ds = xr.open_dataset("ne120np4_nc3000_Co015_Fi001_PF_nullRR_Nsw010_20171011.nc")

    unstructured_to_scrip('test_mesh.nc', ds['lat'].values, ds['lon'].values, debug=True)

    is_valid = validate_scrip_file('test_mesh.nc')
    if not is_valid:
        print("File may not be suitable for ESMF")

    # CMZ note, ESMF is 1-based, so need to subtract one here
    #ESMF_CELL=91
    #test_specific_cell('test_mesh.nc', ESMF_CELL-1)

    #test_grid_coverage('test_mesh.nc')


