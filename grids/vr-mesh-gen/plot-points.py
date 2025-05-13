import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.patches as mpatches
from cartopy.io import shapereader
from shapely.geometry import LineString, MultiLineString

# Load mesh centers
column_names = ['label', 'longitude', 'latitude']
df = pd.read_csv('tc_mesh_centers.csv', names=column_names)

# Set up the plot
plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
# Focus on US east coast
ax.set_extent([-130, -65, 25, 50], crs=ccrs.PlateCarree())

# Add features
ax.add_feature(cfeature.STATES, edgecolor='gray')

# Create a colormap for the mesh centers using tab10 and cycling through colors
num_meshes = len(df)
tab10 = plt.cm.tab10

# Create color assignments by cycling through the tab10 colors
mesh_colors = []
for i in range(num_meshes):
    # This loops through indices 0-9 and then starts over
    # When we have 14 meshes, we'll use colors 0-9 and then 0-3 again
    color_idx = i % 10
    mesh_colors.append(tab10(color_idx/10))

# Create a KD-tree for efficient nearest-neighbor lookups
mesh_points = np.array([[row['longitude'], row['latitude']] for _, row in df.iterrows()])
tree = cKDTree(mesh_points)

# Plot mesh centers
for index, row in df.iterrows():
    plt.plot(row['longitude'], row['latitude'], 'o', color=mesh_colors[index],
             markersize=8, transform=ccrs.PlateCarree())
    plt.text(row['longitude'] + 0.5, row['latitude'], f"{int(row['label']):03d}",
             transform=ccrs.PlateCarree(), fontsize=9)

# Function to check if point is within the US east coast
def is_in_east_coast(lon, lat):
    return -130 <= lon <= -65 and 25 <= lat <= 45

# Get coastline data
coastline_fname = shapereader.natural_earth(category='physical', name='coastline', resolution='10m')
coastline_reader = shapereader.Reader(coastline_fname)

# Process coastline geometries and color segments based on nearest mesh
for record in coastline_reader.records():
    geom = record.geometry

    # Skip geometries outside our area
    if not (geom.bounds[0] <= -65 and geom.bounds[2] >= -85 and
            geom.bounds[1] <= 45 and geom.bounds[3] >= 25):
        continue

    # Process the geometry
    if isinstance(geom, LineString):
        lines = [geom]
    elif isinstance(geom, MultiLineString):
        lines = list(geom.geoms)
    else:
        continue  # Skip other geometry types

    # Process each line
    for line in lines:
        # Get coordinates
        coords = np.array(line.coords)

        # Only process if we have enough coordinates
        if len(coords) < 2:
            continue

        # Create segments and find closest mesh for midpoint of each segment
        for i in range(len(coords) - 1):
            p1 = coords[i]
            p2 = coords[i+1]

            # Calculate midpoint
            mid_x = (p1[0] + p2[0]) / 2
            mid_y = (p1[1] + p2[1]) / 2

            # Only process points in our area of interest
            if not is_in_east_coast(mid_x, mid_y):
                continue

            # Find closest mesh
            _, idx = tree.query([mid_x, mid_y])

            # Plot the segment with appropriate color
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                    color=mesh_colors[idx], linewidth=2,
                    transform=ccrs.PlateCarree(), zorder=10)

# Add a legend
legend_patches = [mpatches.Patch(color=mesh_colors[i],
                                label=f"Mesh {int(row['label']):03d}")
                for i, (_, row) in enumerate(df.iterrows())]
ax.legend(handles=legend_patches, loc='upper left',
         bbox_to_anchor=(1.05, 1), borderaxespad=0.)

plt.title('US East Coast Colored by Closest Weather Model Mesh')
plt.tight_layout()
plt.savefig('colored_coastline_map.png', dpi=300, bbox_inches='tight')