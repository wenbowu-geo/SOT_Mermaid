from scipy.interpolate import RegularGridInterpolator
import numpy as np
from read_etopo1 import read_etopo1_netcdf
import pandas as pd
import argparse

# Set up command-line argument parsing to allow user-defined inputs
parser = argparse.ArgumentParser(description="Extract seafloor depth from ETOPO1 data at specified station locations.")

# Define the directory where ETOPO1 data files are stored
parser.add_argument('--dir_etopo1', type=str, default="/path/to/etopo1/files/", 
                    help="Path to the directory containing ETOPO1 bathymetry data files.")

# Specify the north-south extension range for depth searching (in meters)
parser.add_argument('--extend_m', type=float, default=50000.0, 
                    help="Distance (in meters) to extend the search range north and south from each station \
                    location when determining the deepest seafloor depth.")

args = parser.parse_args()
extend_m = args.extend_m

# Convert extend_m (meters) to degrees of latitude
extend_deg = extend_m / 111320.0  # 1 degree latitude ≈ 111.32 km

# File paths
dir_etopo1 = args.dir_etopo1  
etopo_file = dir_etopo1 + "/etopo5.grd"  

sta_loc_file = "data/loc_Mermaid.txt"
loc_seafloor_cutoff_file = "data/loc_seafloor_cutoff.txt"

# Read the first line (header)
with open(sta_loc_file, "r") as f:
    header = f.readline().strip()  

# Read station data (skipping header row)
data = pd.read_csv(sta_loc_file, delim_whitespace=True, skiprows=1, names=["sta", "lat", "lon", "parking_depth"])

# Convert longitude to 0-360 format if needed
data["lon_360"] = data["lon"] % 360  

# Read ETOPO1 data
lon_etopo1, lat_etopo1, ele_etopo1 = read_etopo1_netcdf(etopo_file, np.min(data["lon_360"].values))

# Create an interpolator for structured grid data
interpolator = RegularGridInterpolator(
    (lat_etopo1, lon_etopo1),  
    ele_etopo1,  
    method="linear",  
    bounds_error=False,  
    fill_value=np.nan  
)

# Initialize depth array
depth_etopo1 = np.zeros(len(data))

# Find maximum depth within extended latitude range
for i in range(len(data)):
    lat = data["lat"].iloc[i]
    lon = data["lon_360"].iloc[i]

    # Define latitude range
    lat_min = lat - extend_deg
    lat_max = lat + extend_deg

    # Get all latitude points within this range (matching the ETOPO1 grid)
    lat_indices = (lat_etopo1 >= lat_min) & (lat_etopo1 <= lat_max)
    lat_grid_points = lat_etopo1[lat_indices]

    # Generate query points for all latitudes at the fixed longitude
    query_points = np.array([[lat_pt, lon] for lat_pt in lat_grid_points])

    # Get depths for all these latitude points
    depths = interpolator(query_points)

    # Store the deepest seafloor depth
    depth_etopo1[i] = np.nanmin(depths)  

# Append "ELVE" to the header
header = header + " ELVE"  
columns = header.split()
col_widths = [max(len(col), 10) for col in columns]  

# Construct format string dynamically
fmt_str = " ".join([f"{{:<{w}}}" for w in col_widths]) + "\n"

# Write output file
with open(loc_seafloor_cutoff_file, "w") as f:
    f.write(fmt_str.format(*columns))  

    for i in range(len(data)):
        f.write(fmt_str.format(
            data["sta"].iloc[i], 
            f"{data['lat'].iloc[i]:.4f}", 
            f"{data['lon'].iloc[i]:.4f}", 
            f"{data['parking_depth'].iloc[i]:.1f}", 
            f"{depth_etopo1[i]:.2f}"
        ))

print(f"Output written to {loc_seafloor_cutoff_file}")

