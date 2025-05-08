import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

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

# Function to subtract some number of hours from date string in format 'YYYYMMDDHH'
def subtract_hours(date_str, hours):
    """
    Subtract `hours` from a timestamp string YYYYMMDDHH
    and return the result in the same format.
    """
    dt = datetime.strptime(date_str, '%Y%m%d%H')
    return (dt - timedelta(hours=hours)).strftime('%Y%m%d%H')

# Define a list of offsets we want, this means "include" times NN offset hours
# before the landfall time in the CSV file
# if empty, just reproduce the landfall time variable
offsets = [12, 24]

# Define the year range for filtering landfall events
min_year = 1988
max_year = 2023

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

    # Get the list of dates for this mesh center
    current_dates = filtered_data['Ttime'].tolist()

    # Build a list of dates associated with this landfall and write to file
    if offsets:
        expanded_dates = []
        for date in current_dates:
            # loop offsets largest→smallest
            for h in sorted(offsets, reverse=True):
                expanded_dates.append(subtract_hours(date, h))
            # then include the original date
            expanded_dates.append(date)

        # write expanded dates to file
        with open(filename, 'w') as f:
            for d in expanded_dates:
                f.write(f"{d}\n")

        all_dates.extend(expanded_dates)
    else:
        # if offsets is empty, do nothing special
        filtered_data.to_csv(filename, index=False, header=False)
        all_dates.extend(current_dates)

# We also need to create a single consolidated file containing all unique landfall dates
# This is used for land spinup or single grid runs
# First, remove duplicates (when landfalls occurred on same date in different regions)
# Then sort chronologically
unique_sorted_dates = sorted(set(all_dates))

# Write the consolidated dates to a single file (useful for CLM spinup)
with open('all_dates_sorted.txt', 'w') as f:
    for date in unique_sorted_dates:
        f.write(f"{date}\n")