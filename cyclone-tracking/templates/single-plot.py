import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
from io import StringIO

# The data as a string
data_str = """
    1718794 289.192922      35.332488       9.367170e+04    5.162286e+01    2075    9       4       0
    1719296 288.796697      35.646568       9.363564e+04    5.226929e+01    2075    9       4       3
    1734055 288.666962      35.890092       9.329712e+04    5.460104e+01    2075    9       4       6
    1734687 288.077623      36.380170       9.324256e+04    5.326963e+01    2075    9       4       9
    1733838 287.750480      36.404635       9.289647e+04    5.132148e+01    2075    9       4       12
    1733977 287.560605      36.559479       9.320971e+04    4.795106e+01    2075    9       4       15
    1749539 287.644538      36.778779       9.329451e+04    4.842957e+01    2075    9       4       18
    1750158 287.317569      37.285464       9.332232e+04    4.924231e+01    2075    9       4       21
    1766517 287.236480      37.834969       9.349752e+04    4.893792e+01    2075    9       5       0
    1781940 286.789770      38.384334       9.385073e+04    5.221305e+01    2075    9       5       3
    1798010 286.364341      38.954087       9.414823e+04    5.214890e+01    2075    9       5       6
    1813076 285.758295      39.696268       9.457373e+04    4.647403e+01    2075    9       5       9
    1829643 285.202268      40.238287       9.541409e+04    3.808104e+01    2075    9       5       12
    1844918 284.491297      40.818280       9.673873e+04    2.978595e+01    2075    9       5       15
    1844596 283.975961      41.333206       9.779727e+04    2.163022e+01    2075    9       5       18
    1860174 283.807899      41.811990       9.877789e+04    1.650366e+01    2075    9       5       21
    1876480 283.755865      42.293855       9.926745e+04    1.415181e+01    2075    9       6       0
    483255  285.413148      44.943325       9.941458e+04    2.032780e+01    2075    9       6       3
"""

# Read the data into a Pandas DataFrame
data = pd.read_csv(StringIO(data_str), sep="\s+", header=None, names=["ID", "Longitude", "Latitude","psl","wind","yyyy","mm","dd","hh"])

# Extract longitude and latitude for plotting
longitude = data['Longitude']
latitude = data['Latitude']

# Create the map plot
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAKES, alpha=0.5)
ax.add_feature(cfeature.RIVERS)

# Plot the cyclone path
plt.plot(longitude, latitude, marker='o', linestyle='-', color='b', transform=ccrs.PlateCarree())

# Add gridlines and labels
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

# Zoom the map to the region of interest
ax.set_extent([min(longitude) - 5, max(longitude) + 5, min(latitude) - 5, max(latitude) + 5], crs=ccrs.PlateCarree())

plt.show()
