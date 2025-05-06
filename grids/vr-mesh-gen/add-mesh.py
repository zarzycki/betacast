import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Function to calculate the Haversine distance between two lat/lon pairs
def haversine(lon1, lat1, lon2, lat2):
    # degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371.
    return c * r

# Define column names and load landfall file using pandas
# Ttime is the "time of initialization" which is N days before actual LFtime (defined by LF code)
# That is the time the model should be initialized in hindcast mode
column_names = ['stormID', 'LFtime', 'LFlon','LFlat','LFpres','LFwind','Tlon','Tlat','Ttime','LFbasin']
df = pd.read_csv('LF.ibtracs-1980-2019-GLOB.v4.txt', names=column_names)

# Only keep regions of landfall corresponding to USA
filtered_df = df[df['LFbasin'] > 0].copy()
print(filtered_df)
print(f"Number of rows kept: {filtered_df.shape[0]}")

# Define columns and load TC mesh centers information
column_names2 = ['label', 'longitude', 'latitude']
df2 = pd.read_csv('tc_mesh_centers.csv', names=column_names2)
# Convert 'label' to integer and then format it to be zero-padded to three digits
df2['label'] = df2['label'].astype(int).apply(lambda x: f"{x:03}")

# Initialize a new column in the TC landfalls for the labels
filtered_df['closest_label'] = ''

# Iterate over each row in filtered_df to calculate how far each center point from df2 is
# Find the minimum and assign it
for index, row in filtered_df.iterrows():
    min_distance = np.inf
    closest_label = None
    for _, row2 in df2.iterrows():
        distance = haversine(row['LFlon'], row['LFlat'], row2['longitude'], row2['latitude'])
        if distance < min_distance:
            min_distance = distance
            closest_label = row2['label']
    filtered_df.at[index, 'closest_label'] = closest_label

# Histogram sorted by labels alphabetically
label_counts = filtered_df['closest_label'].value_counts().sort_index()
label_counts.plot(kind='bar')
plt.xlabel('Label')
plt.ylabel('Frequency')
plt.title('Histogram of Closest Label')
plt.xticks(rotation=0)  # Ensures the x-axis labels are not rotated and are easily readable
plt.savefig('histogram.png', dpi=300, bbox_inches='tight')

# Define the year range for filtering landfall events
min_year = 1988
max_year = 2023

# Specify which column contains the timestamp information for output
output_columns = ['Ttime']

# Initialize a list to collect all dates across all mesh centers
all_dates = []

# Process each unique mesh center label (previously assigned in the distance calculation)
for label in filtered_df['closest_label'].unique():

    # Filter data to only rows associated with the current mesh center label
    # ex: get all rows with "001"
    filtered_data = filtered_df[filtered_df['closest_label'] == label][output_columns]

    # Convert Ttime to string format
    # Filter events by year, assuming Ttime format starts with YYYY
    # This keeps only events that occurred within the specified year range
    filtered_data['Ttime'] = filtered_data['Ttime'].astype(str)
    filtered_data = filtered_data[filtered_data['Ttime'].apply(lambda x: min_year <= int(x[:4]) <= max_year)]

    # Generate filename for this mesh center label's date output file
    filename = f"dates.index.{label}.txt"

    # Write the filtered dates to a separate file for each mesh center
    # These individual files can be used for region-specific inits
    filtered_data.to_csv(filename, index=False, header=False)

    # Collect all dates
    all_dates.extend(filtered_data['Ttime'].tolist())

# Create a single consolidated file containing all unique landfall dates
# First, remove duplicates (when landfalls occurred on same date in different regions)
# Then sort chronologically
unique_sorted_dates = sorted(set(all_dates))

# Write the consolidated dates to a single file (useful for CLM spinup)
with open('all_dates_sorted.txt', 'w') as f:
    for date in unique_sorted_dates:
        f.write(f"{date}\n")