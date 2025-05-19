#!/usr/bin/env python3

import math
import gzip

kts_to_ms = 0.514444
nm_to_km = 1.852

def parse_EBT(line):
    """
    Parse a single Extended Best Track format line and return structured data
    Uses column-based parsing for fixed-width fields, then space-based for remainder
    """
    if not line.strip():
        return None

    # Print raw line for debugging
    print(f"Raw EBT line: {repr(line)}")
    print(f"Line length: {len(line)} characters")

    # Remove newline but preserve the full line
    line = line.rstrip('\n')

    # Define column positions based on EBT format
    # https://rammb2.cira.colostate.edu/wp-content/uploads/2020/11/EBTRK_README_v3.0.2_16-Oct-2022.txt
    storm_id     = line[0:9].strip()      # AL132019
    storm_name   = line[9:21].strip()     # LORENZO
    date_time    = line[21:27].strip()    # MMDDHH
    year         = line[27:32].strip()
    lat          = line[32:39].strip()
    lon          = line[39:46].strip()
    vmax         = line[46:49].strip()    # max wind (kt)
    mslp         = line[49:54].strip()    # min pressure (hPa)
    rmw          = line[54:58].strip()    # radius of max wind (nm)
    eye          = line[58:62].strip()    # eye radius (nm)
    outer_p      = line[62:67].strip()    # outer closed isobar pressure (hPa)
    r8           = line[67:71].strip()    # outer closed isobar radius (nm)
    r34          = line[71:84].strip()    # 12 chat 34kt quadrants (nm) [NE, SE, SW, NW]
    r50          = line[84:97].strip()    # 12 chat 50kt quadrants (nm) [NE, SE, SW, NW]
    r64          = line[97:110].strip()   # 12 chat 64kt quadrants (nm) [NE, SE, SW, NW]
    storm_type   = line[111:112].strip()  # e.g., 'L'

    print("Column-parsed fields:")
    print(f"  Storm ID:      '{storm_id}'")
    print(f"  Storm Name:    '{storm_name}'")
    print(f"  Date/Time:     '{date_time}'")
    print(f"  Year:          '{year}'")
    print(f"  Latitude:      '{lat}'")
    print(f"  Longitude:     '{lon}'")
    print(f"  Max Wind:      '{vmax}' kt")
    print(f"  MSLP:          '{mslp}' hPa")
    print(f"  RMW:           '{rmw}' nm")
    print(f"  Eye rad:       '{eye}' nm")
    print(f"  Outer P:       '{outer_p}' hPa")
    print(f"  R8:            '{r8}' nm")
    print(f"  r34:           '{r34}'")
    print(f"  r50:           '{r50}'")
    print(f"  r64:           '{r64}'")
    print(f"  Storm Type:    '{storm_type}'")

    # Function for turning 12 char quadrant strings into a list x 4
    def parse_quadrants(s):
        s = s.rjust(12)
        return [float(val) if val.strip() not in ['', '-99'] else None
                for val in [s[i:i+3] for i in range(0, 12, 3)]]

    r34_quads = parse_quadrants(r34)
    r50_quads = parse_quadrants(r50)
    r64_quads = parse_quadrants(r64)

    # Parse month, day, hour from date_time
    if len(date_time) == 6:
        month = date_time[:2]
        day = date_time[2:4]
        hour = date_time[4:6]
    else:
        print(f"WARNING: Unexpected date_time format: {date_time}")
        month = day = hour = "00"

    # Extract basin and storm number from storm_id
    if len(storm_id) >= 4:
        basin_code = storm_id[:2]  # AL, EP, CP, WP
        storm_num = storm_id[2:4]  # 01-99
    else:
        print(f"WARNING: Unexpected storm_id format: {storm_id}")
        basin_code = "AL"
        storm_num = "01"

    parsed_data = {
        'storm_id': storm_id,
        'storm_name': storm_name,
        'year': year,
        'month': month,
        'day': day,
        'hour': hour,
        'lat': float(lat),
        'lon': float(lon),
        'max_wind_kt': float(vmax),
        'min_pressure': int(mslp),
        'rmw_nm': float(rmw),
        'eye_diam': float(eye),
        'outer_pressure': float(outer_p),
        'outer_radius_nm': float(r8),
        'wind_34': r34_quads,
        'wind_50': r50_quads,
        'wind_64': r64_quads,
        'storm_type_code': storm_type,
        'basin_code': basin_code,
        'storm_num': storm_num
    }

    print(f"Parsed EBT data: {parsed_data}")
    print("-" * 80)

    return parsed_data


def convert_to_tcvitals(ebt_data, prev_lat=None, prev_lon=None):
    """
    Convert parsed EBT data to TCVitals format
    """
    print(f"\nConverting to TCVitals: {ebt_data['storm_name']} {ebt_data['year']}{ebt_data['month']}{ebt_data['day']} {ebt_data['hour']}00")

    # Calculate motion (for now, just set to missing values)
    motion_dir, motion_speed = -99, -99

    # Convert basin code to single letter
    basin_letter_map = {'AL': 'L', 'EP': 'E', 'CP': 'C', 'WP': 'W'}
    basin_letter = basin_letter_map.get(ebt_data['basin_code'], 'L')

    # Format storm identifier (e.g., 04L)
    storm_identifier = f"{ebt_data['storm_num']}{basin_letter}"

    # Convert wind from kt to m/s
    max_wind_ms = int(ebt_data['max_wind_kt'] * kts_to_ms)

    # Format coordinates
    lat_tenths = f"{int(ebt_data['lat'] * 10):03d}"
    lat_flag = 'N'
    lon_tenths = f"{int(ebt_data['lon'] * 10):04d}"
    lon_flag = 'W'

    # Convert radii from nm to km and format
    def convert_radii(radii_nm):
        converted = []
        for r in radii_nm:
            if r is None or r == -99 or r == 0:
                converted.append('-999')
            else:
                # Convert nm to km (multiply by 1.852)
                km = int(r * 1.852)
                converted.append(f"{km:04d}")
        return converted

    wind_34_km = convert_radii(ebt_data['wind_34'])
    wind_50_km = convert_radii(ebt_data['wind_50'])
    wind_64_km = convert_radii(ebt_data['wind_64'])

    print(f"Wind radii converted to km:")
    print(f"  34kt: {wind_34_km}")
    print(f"  50kt: {wind_50_km}")
    print(f"  64kt: {wind_64_km}")

    # Convert other measurements
    if ebt_data['outer_radius_nm'] == '-99' or float(ebt_data['outer_radius_nm']) < 0:
        outer_radius_str = '-999'
    else:
        outer_radius_km = int(float(ebt_data['outer_radius_nm']) * nm_to_km)
        outer_radius_str = f"{outer_radius_km:04d}"

    if ebt_data['rmw_nm'] == '-99' or float(ebt_data['rmw_nm']) < 0:
        rmw_str = '-99'
    else:
        rmw_km = int(float(ebt_data['rmw_nm']) * nm_to_km)
        rmw_str = f"{rmw_km:03d}"

    if ebt_data['min_pressure'] < 0:
        min_pressure_str = "-999"
    else:
        min_pressure_str = f"{ebt_data['min_pressure']:04d}"

    if ebt_data['outer_pressure'] < 0:
        env_pressure = "-999"
    else:
        env_pressure = f"{int(ebt_data['outer_pressure']):04d}"

    # Build TCVitals line with exact spacing
    tcvital_line = (
        f"NHC  {storm_identifier} {ebt_data['storm_name']:<9} "
        f"{ebt_data['year']}{ebt_data['month']}{ebt_data['day']} {ebt_data['hour']}00 "
        f"{lat_tenths}{lat_flag} {lon_tenths}{lon_flag} "
        f"{-99:03d} {-99:03d} "
        f"{min_pressure_str} {env_pressure} "
        f"{outer_radius_str} {max_wind_ms:02d} "
        f"{rmw_str} "
        f"{wind_34_km[0]} {wind_34_km[1]} {wind_34_km[2]} {wind_34_km[3]} "
        f"X "   # Vortex depth, X for now
        f"{wind_50_km[0]} {wind_50_km[1]} {wind_50_km[2]} {wind_50_km[3]} "
        f"-9 -99N -999W "  # Default forecast values
        f"{wind_64_km[0]} {wind_64_km[1]} {wind_64_km[2]} {wind_64_km[3]} "
    )

    print(f"Generated TCVitals line:")
    print(f"  {tcvital_line}")
    print(f"  Length: {len(tcvital_line)} characters")

    return tcvital_line

def convert_ebtrk_to_tcvitals(input_file, filter_year=None):
    """
    Convert Extended Best Track format to TCVitals format

    Args:
        input_file: Path to the EBT input file
        filter_year: Optional year to filter for (e.g., 2002). If None, processes all years.
    """
    print(f"Converting {input_file} to TCVitals format")
    if filter_year:
        print(f"Filtering for year: {filter_year}")
    print("=" * 80)

    # Read the input file
    try:
        if input_file.endswith('.gz'):
            print("Reading a gzip file!")
            with gzip.open(input_file, 'rt') as f:
                lines = f.readlines()
        else:
            print("Reading an ASCII/text file!")
            with open(input_file, 'r') as f:
                lines = f.readlines()
    except FileNotFoundError:
        print(f"ERROR: Could not find input file: {input_file}")
        return []

    print(f"Read {len(lines)} lines from input file")

    tcvitals_lines = []
    prev_lat, prev_lon = None, None
    total_processed = 0
    filtered_count = 0

    for i, line in enumerate(lines):
        # Parse the EBT line first
        ebt_data = parse_EBT(line)
        if ebt_data is None:
            print(f"Skipping invalid line {i+1}")
            continue

        # Apply year filter if specified
        if filter_year is not None:
            if int(ebt_data['year']) != filter_year:
                filtered_count += 1
                continue

        print(f"\n--- Processing line {i+1} (Year: {ebt_data['year']}) ---")

        # Convert to TCVitals
        tcvital_line = convert_to_tcvitals(ebt_data, prev_lat, prev_lon)
        tcvitals_lines.append(tcvital_line)

        # Update previous position for motion calculation
        prev_lat, prev_lon = ebt_data['lat'], ebt_data['lon']
        total_processed += 1

    if filter_year:
        print(f"\nFiltering results: {filtered_count} records filtered out, {total_processed} records kept for year {filter_year}")

    return tcvitals_lines


if __name__ == "__main__":

    # Filter for year (set to None for full dataset)
    target_year = None

    # EBT file
    input_filename = 'EBTRK_AL_final_1851-2021_new_format_02-Sep-2022-1.txt.gz'

    if target_year:
        output_filename = f'vitals_{target_year}.txt'
    else:
        output_filename = 'vitals_ALL.txt'

    print("Hurricane Track Converter - Extended Best Track to TCVitals")
    print("=" * 60)

    tcvitals_data = convert_ebtrk_to_tcvitals(input_filename, filter_year=target_year)

    if tcvitals_data:
        print(f"\n\nFINAL RESULTS:")
        print("=" * 60)
        print("TCVitals format output:")
        print("-" * 155)
        for line in tcvitals_data:
            print(line)

        # Save to file
        with open(output_filename, 'w') as f:
            for line in tcvitals_data:
                f.write(line + '\n')

        print(f"\nConverted {len(tcvitals_data)} records for year {target_year}")
        print(f"Output saved to '{output_filename}'")
    else:
        print(f"No valid records were converted for year {target_year}")