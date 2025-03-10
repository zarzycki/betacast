import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import pandas as pd

column_names = ['label', 'longitude', 'latitude']
df = pd.read_csv('tc_mesh_centers.csv', names=column_names)

plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Zoom in on US
ax.set_extent([-130, -65, 20, 50], crs=ccrs.PlateCarree())

ax.coastlines()
ax.add_feature(cfeature.STATES)

for index, row in df.iterrows():
    plt.plot(row['longitude'], row['latitude'], 'bo', transform=ccrs.Geodetic())
    plt.text(row['longitude'] + 0.5, row['latitude'], f"{int(row['label']):03d}", transform=ccrs.Geodetic())

plt.savefig('map_plot.png', dpi=300, bbox_inches='tight')