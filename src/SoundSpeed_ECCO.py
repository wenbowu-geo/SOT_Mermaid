import argparse
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy import arange, exp
import xarray as xr
import gsw
import great_circle_calculator.great_circle_calculator as gcc
import itertools
import warnings
from datetime import datetime, timedelta, timezone


warnings.filterwarnings('ignore')


# --- Utility Functions ---

def filledges(a):
    """Fill in coastal and bottom points, ignoring edge points."""
    nx, ny, nz = a.shape
    b = np.copy(a)
    for i, j, k in itertools.product(range(1, nx-1), range(1, ny-1), range(nz)):
        if np.isnan(a[i, j, k]) and np.any(~np.isnan(a[i-1:i+2, j-1:j+2, k])):
            b[i, j, k] = np.nanmean(a[i-1:i+2, j-1:j+2, k])
    for i, j, k in itertools.product(range(nx), range(ny), range(nz-1)):
        if ~np.isnan(a[i, j, k]) and np.isnan(a[i, j, k+1]):
            b[i, j, k+1] = a[i, j, k]
    return b

def calculate_mean_st(year0_mean, year1_mean, dir_ECCO):
    """
    Calculate mean salinity and temperature over a range of years.
    Parameters:
        year0_mean (int): Start year for calculating mean.
        year1_mean (int): End year for calculating mean.
        dir_ECCO (str): Directory path where the ECCO data files are located.
        meanST_exist (bool): Flag to check if pre-computed mean files exist. Default is False.
    Returns:
        SALT_mean (np.array): Mean salinity values.
        THETA_mean (np.array): Mean temperature values.
    """
    # Open the first month's dataset to initialize arrays
    f_S = xr.open_mfdataset(dir_ECCO + 'SALT_' + str(year0_mean) + '_01.nc')
    f_T = xr.open_mfdataset(dir_ECCO + 'THETA_' + str(year0_mean) + '_01.nc')
    # Initialize mean arrays with the shape of the SALT and THETA data
    SALT_mean = np.zeros(f_S.SALT.values.shape)
    THETA_mean = np.zeros(f_T.THETA.values.shape)
    # Define number of years and months
    nyears_ECCO = year1_mean - year0_mean + 1
    years_ECCO = np.arange(year0_mean, year1_mean + 1)
    print(years_ECCO)
    # Loop through the years and months to accumulate values for mean calculation
    for iyear in years_ECCO:
        for imonth in range(1, 13):
            file_SALT_this_month = dir_ECCO + 'SALT_' + str(iyear) + '_' + str(imonth).zfill(2) + '.nc'
            file_THETA_this_month = dir_ECCO + 'THETA_' + str(iyear) + '_' + str(imonth).zfill(2) + '.nc'
            # Load SALT and THETA data for the current month
            fl_tmp = xr.open_mfdataset(file_SALT_this_month)
            SALT_mean += fl_tmp.SALT.values
            fl_tmp = xr.open_mfdataset(file_THETA_this_month)
            THETA_mean += fl_tmp.THETA.values
    # Calculate the mean values over all months and years
    SALT_mean /= (nyears_ECCO * 12.0)
    THETA_mean /= (nyears_ECCO * 12.0)
    # Mask the missing data (zero values)
    SALT_mean[SALT_mean == 0.0] = np.nan
    THETA_mean[THETA_mean == 0.0] = np.nan
    # Assign the computed mean values to the dataset objects
    f_S.SALT.values = SALT_mean
    f_T.THETA.values = THETA_mean
    print(f_S.SALT.values[0, 30, 100, 100], SALT_mean[0, 30, 100, 100])
    # Save the results to netCDF files
    f_S.to_netcdf('SALT.ECCOmean.nc')
    f_T.to_netcdf('THETA.ECCOmean.nc')
    # Apply edge filling for the data (for spatial continuity)
    a_tmp = np.empty(SALT_mean[0, :, :, :].shape)
    a_tmp = SALT_mean[0, :, :, :]
    SALT_extended = filledges(a_tmp)
    a_tmp = THETA_mean[0, :, :, :]
    THETA_extended = filledges(a_tmp)
    # Get the latitude, longitude, and depth information
    lat_ECCO = fl_tmp.latitude.values
    lon_ECCO = fl_tmp.longitude.values
    dep_ECCO = -fl_tmp.Z.values
    # Assign the extended data to the dataset
    f_S.SALT.values[0, :, :, :] = SALT_extended
    f_T.THETA.values[0, :, :, :] = THETA_extended
    print(f_S.SALT.values[0, 30, 100, 100], SALT_mean[0, 30, 100, 100])
    # Save the extended results to netCDF files
    f_S.to_netcdf('SALT.ECCOmean.extended.nc')
    f_T.to_netcdf('THETA.ECCOmean.extended.nc')
    # Close the datasets
    f_S.close()
    f_T.close()
    return SALT_mean, THETA_mean,lon_ECCO,lat_ECCO,dep_ECCO

def fillbtm2(a):
    """Fill missing bottom points with nearest valid values."""
    nx, nz = a.shape
    b = np.copy(a)
    for i in range(nx):
        for k in range(nz-1, 0, -1):
            if ~np.isnan(b[i, k]) and np.isnan(b[i, k-1]):
                b[i, k-1] = b[i, k]
    return b


def ip(frac, p1_extend, p2_extend):
    """Calculate intermediate points along a great-circle path."""
    return gcc.intermediate_point(p1_extend, p2_extend, frac)


def c(T, S, lon, lat, z):
    """Calculate sound speed from in-situ temperature and salinity."""
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    insitu_t = gsw.t_from_CT(SA, T, p)
    return gsw.sound_speed_t_exact(SA, insitu_t, p)


def rho(T, S, lon, lat, z):
    """Calculate in-situ density from salinity and temperature."""
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    return gsw.rho_t_exact(SA, T, p)

def load_mean_st(dir_ECCO):
    """Load precomputed mean salinity and temperature."""
    f_S = xr.open_mfdataset(f"{dir_ECCO}/SALT.ECCOmean.nc")
    f_T = xr.open_mfdataset(f"{dir_ECCO}/THETA.ECCOmean.nc")
    SALT_extended = f_S.SALT.values[0, :, :, :]
    THETA_extended = f_T.THETA.values[0, :, :, :]
    lat_ECCO = f_S.latitude.values
    lon_ECCO = f_S.longitude.values
    dep_ECCO = -f_S.Z.values  # Depths are negative in ECCO
    return SALT_extended, THETA_extended, lat_ECCO, lon_ECCO, dep_ECCO

def compute_great_circle_path(psrc, prec, extend_m, path_dx_m):
    """
    Compute points along the great-circle path and calculate signed distances relative to the receiver.
    Parameters:
        psrc (tuple): Longitude and latitude of the source (lon_src, lat_src).
        prec (tuple): Longitude and latitude of the receiver (lon_rec, lat_rec).
        extend_m (float): Distance to extend the path beyond source and receiver, in meters.
        path_dx_m (float): Spacing between points along the path, in meters.
    Returns:
        d0 (float): Distance between the source and receiver, in meters.
        grid_path (array): Array of latitude and longitude along the path [(lat1, lon1), (lat2, lon2), ...].
        lat_path (array): Latitudes along the path.
        lon_path (array): Longitudes along the path.
        x_path (array): Signed distances relative to the receiver (positive toward source, negative away).
    """
    # Compute source-receiver distance and bearings
    d0 = gcc.distance_between_points(psrc, prec, unit='meters')
    crs1, crs2 = gcc.bearing_at_p1(psrc, prec), gcc.bearing_at_p1(prec, psrc)
    # Extend path beyond source and receiver
    p1_extend = gcc.point_given_start_and_bearing(psrc, crs2, extend_m, unit='meters')
    p2_extend = gcc.point_given_start_and_bearing(prec, crs1, extend_m, unit='meters')
    # Interpolate points along the great-circle path
    path_xnpts = int((d0 + 2 * extend_m) / path_dx_m) + 1
    fr = np.linspace(0, 1, path_xnpts)
    grid_path = np.array([ip(f, p2_extend, p1_extend) for f in fr])
    lat_path, lon_path = grid_path[:, 1], grid_path[:, 0]
    # Compute signed distances relative to the receiver
    x_path = np.zeros(len(grid_path))
    for i, point in enumerate(grid_path):
        distance = gcc.distance_between_points(point, prec, unit='meters')
        x_path[i] = distance if point[1] > prec[1] else -distance  # Positive toward source
    return d0, grid_path, lat_path, lon_path, x_path,p1_extend,p2_extend

def interpolate_properties_along_path(lon_cut_ECCO, lat_cut_ECCO, c_cut_ECCO, rho_cut_ECCO, 
                                      T_cut_ECCO, S_cut_ECCO, grid_path, dep_ECCO, 
                                      fillbtm2, path_depth0_m, path_depth1_m, path_ddepth_m):
    """
    Interpolates properties (c, rho, theta, salt) along the selected path at given depths.
    Parameters:
    lon_cut_ECCO (array): Longitudes of ECCO data points.
    lat_cut_ECCO (array): Latitudes of ECCO data points.
    c_cut_ECCO (array): Sound speed data.
    rho_cut_ECCO (array): Density data.
    T_cut_ECCO (array): Temperature data.
    S_cut_ECCO (array): Salinity data.
    grid_path (array): Path coordinates for interpolation.
    dep_ECCO (array): Depth values corresponding to ECCO data.
    fillbtm2 (function): Function to fill in missing values at the bottom.
    path_depth0_m (float): Starting depth.
    path_depth1_m (float): Ending depth.
    path_ddepth_m (float): Depth interval.
    Returns:
    c_path (array): Interpolated sound speed along the path at each depth.
    rho_path (array): Interpolated density along the path at each depth.
    theta_path (array): Interpolated temperature along the path at each depth.
    salt_path (array): Interpolated salinity along the path at each depth.
    """
    # Initialize arrays to store the results (horizontal points x depth points)
    path_xnpts = len(grid_path)
    c_pathx_ECCOdep = np.full((path_xnpts, dep_ECCO.shape[0]), np.nan)
    rho_pathx_ECCOdep = np.full((path_xnpts, dep_ECCO.shape[0]), np.nan)
    theta_pathx_ECCOdep = np.full((path_xnpts, dep_ECCO.shape[0]), np.nan)
    salt_pathx_ECCOdep = np.full((path_xnpts, dep_ECCO.shape[0]), np.nan)
    # Horizontal interpolation for each depth
    for k in range(dep_ECCO.shape[0]):
        itc = RGI((lon_cut_ECCO, lat_cut_ECCO), c_cut_ECCO[:,:,k])
        itr = RGI((lon_cut_ECCO, lat_cut_ECCO), rho_cut_ECCO[:,:,k])
        itt = RGI((lon_cut_ECCO, lat_cut_ECCO), T_cut_ECCO[:,:,k])
        its = RGI((lon_cut_ECCO, lat_cut_ECCO), S_cut_ECCO[:,:,k])
        c_pathx_ECCOdep[:,k] = itc(grid_path)
        rho_pathx_ECCOdep[:,k] = itr(grid_path)
        theta_pathx_ECCOdep[:,k] = itt(grid_path)
        salt_pathx_ECCOdep[:,k] = its(grid_path)
    # Define depth path for interpolation
    depth_path = np.arange(path_depth0_m, path_depth1_m, path_ddepth_m)
    ndepth_path = depth_path.shape[0]
    # Fill in missing bottom values
    c_fillgap = fillbtm2(c_pathx_ECCOdep)
    rho_fillgap = fillbtm2(rho_pathx_ECCOdep)
    theta_fillgap = fillbtm2(theta_pathx_ECCOdep)
    salt_fillgap = fillbtm2(salt_pathx_ECCOdep)
    # Initialize arrays to store interpolated results for each depth
    c_path = np.full((path_xnpts, ndepth_path), np.nan)
    rho_path = np.full((path_xnpts, ndepth_path), np.nan)
    theta_path = np.full((path_xnpts, ndepth_path), np.nan)
    salt_path = np.full((path_xnpts, ndepth_path), np.nan)
    # Depth-wise interpolation for each property
    for i in range(path_xnpts):
        itc = interp1d(dep_ECCO, c_fillgap[i, :], fill_value="extrapolate")
        c_path[i, :] = itc(depth_path)
        itr = interp1d(dep_ECCO, rho_fillgap[i, :], fill_value="extrapolate")
        rho_path[i, :] = itr(depth_path)
        itt = interp1d(dep_ECCO, theta_fillgap[i, :], fill_value="extrapolate")
        theta_path[i, :] = itt(depth_path)
        its = interp1d(dep_ECCO, salt_fillgap[i, :], fill_value="extrapolate")
        salt_path[i, :] = its(depth_path)
    return depth_path,c_path, rho_path, theta_path, salt_path


def save_results(output_dir,tomo_file_name, grid_file_name, ave1D_file_name,grid_path, x_path, depth_path, 
                 c_path, rho_path, salt_path, theta_path, path_dx_m, path_ddepth_m):
    """
    Save computed results to a file, including average properties over the path.
    Parameters:
        tomo_file_name (str): Name of the output file.
        x_path (array): The x positions along the path (in meters).
        depth_path (array): The depth values along the path (in meters).
        c_path (array): Sound speed along the path and depth.
        rho_path (array): Density along the path and depth.
        salt_path (array): Salinity along the path and depth.
        theta_path (array): Temperature along the path and depth.
        path_dx_m (float): Spacing between points along the horizontal path (in meters).
        path_ddepth_m (float): Spacing between points along the depth path (in meters).
    """
    # Compute averages over distance (axis=0) while ignoring NaNs
    speed_average = np.nanmean(c_path, axis=0)  # Average over distance for each depth
    density_average = np.nanmean(rho_path, axis=0)
    salt_average = np.nanmean(salt_path, axis=0)
    theta_average = np.nanmean(theta_path, axis=0)
    # Write results to file
    head1 = f'{x_path[0]:.2f} {depth_path[0]:.2f} {x_path[-1]:.2f} {depth_path[-1]:.2f} #x_min(m) z_min(m) x_max(m) z_max(m)\n'
    head2 = f'{path_dx_m:.4f} {path_ddepth_m:.4f}                #dx(m) dz(m)\n'
    head3 = f'{len(x_path)} {len(depth_path)}                           #nx nz\n'
    with open(output_dir+'/'+tomo_file_name, "w") as file:
        file.write(head1)
        file.write(head2)
        file.write(head3)
        for i, k in itertools.product(range(len(x_path)), range(len(depth_path))):
            file.write(f'{x_path[i]:.2f} {abs(depth_path[k]):.2f} {c_path[i, k]:.4f}  {rho_path[i, k]:.4f} {salt_path[i, k]:.4f} {theta_path[i, k]:.4f}\n')
    print(f"Results saved to {tomo_file_name}")
    # Save average results to the ave1D file
    with open(output_dir+'/'+ave1D_file_name, "w") as file:
        file.write("# Averages over distance (i) for each depth (k):\n")
        file.write("# Depth(m) Speed(m/s) Density(kg/m^3) Salt PSU Theta(degC)\n")
        for k in range(len(depth_path)):
            file.write(f'{abs(depth_path[k]):.2f} {speed_average[k]:.4f} {density_average[k]:.4f} '
                       f'{salt_average[k]:.4f} {theta_average[k]:.4f}\n')
    print(f"Average results saved to {ave1D_file_name}")
    np.savez(output_dir+'/'+grid_file_name,
             sound_speed=c_path,
             density=rho_path,
             temperature=theta_path,
             salinity=salt_path,
             grid_path=grid_path,
             depth_path=depth_path,
             distance_path=x_path)
    return speed_average, density_average, salt_average, theta_average


# --- Main Script ---
tomo_file_name="tomo_acoustic.txt"
ave1D_file_name="average1D_VpRhoTS.txt"
grid_file_name="property_path.npz"

parser = argparse.ArgumentParser(description="Process source-receiver parameters and model coverage settings.")
parser.add_argument('--psrc', type=float, nargs=2, required=True, help="Source coordinates (longitude, latitude).")
parser.add_argument('--prec', type=float, nargs=2, required=True, help="Receiver coordinates (longitude, latitude).")
parser.add_argument('--extend_m', type=float, default=50000.0, help="Extend source-receiver path at both ends (meters).")
parser.add_argument('--path_dx_m', type=float, default=5000.0, help="Path spacing in x-direction (meters).")
parser.add_argument('--path_depth0_m', type=float, default=0.0, help="Starting depth of the path (meters).")
parser.add_argument('--path_depth1_m', type=float, default=8000.0, help="Ending depth of the path (meters).")
parser.add_argument('--path_ddepth_m', type=float, default=10.0, help="Path spacing in depth (meters).")
parser.add_argument('--meanST_exist', type=bool, default=True, help="Whether mean sound speed exists (True/False).")
parser.add_argument('--year0_mean', type=int, default=1992, help="Starting year for mean calculation.")
parser.add_argument('--year1_mean', type=int, default=1992, help="Ending year for mean calculation.")
parser.add_argument('--dir_ECCO', type=str, default="/path/to/ECCO/files/", help="Directory for ECCO files.")
parser.add_argument('--output_dir', type=str, default="./output", help="Directory to save results.")

args = parser.parse_args()

# Unpack arguments
psrc, prec = tuple(args.psrc), tuple(args.prec)
extend_m, path_dx_m = args.extend_m, args.path_dx_m
path_depth0_m, path_depth1_m = args.path_depth0_m, args.path_depth1_m
path_ddepth_m = args.path_ddepth_m
meanST_exist, year0_mean, year1_mean = args.meanST_exist, args.year0_mean, args.year1_mean
dir_ECCO, output_dir = args.dir_ECCO, args.output_dir

#dir_ECCO="/proj/mazu/wenbowu/GeophyDatabase/Oceanography/ECCOv4_Release4/data/"
#dir_ECCO="~/work/GeophyDatabase/Oceanography/data/"
#dir_ECCO="/Users/wenbowu/work/T_wave/codes/debug/MakeSoundSpeed/"
#year0_mean=1992
#year1_mean=1992
#source and receiver locations
#psrc,prec=(lon,lat),(lon,lat)
#psrc,prec=(-78.949699,-33.446701), (-78.846199,-33.833801)
#extend the source-receiver path at both ends for ensure complete model coverage.
#extend_m=50000.0
#path_dx_m=5000.0
#path_depth0_m=0.0
#path_depth1_m=8000.0
#path_ddepth_m=10.0
#meanST_exist=True
#output_dir="./"

# Load ECCO data
print("Read ECCO data...")
if not meanST_exist:
    SALT_mean, THETA_mean, lon_ECCO,lat_ECCO,dep_ECCO = calculate_mean_st(year0_mean, year1_mean, dir_ECCO)
    SALT_extended = filledges(SALT_mean[0, :, :, :])
    THETA_extended = filledges(THETA_mean[0, :, :, :])
else:
    SALT_extended, THETA_extended, lat_ECCO, lon_ECCO, dep_ECCO = load_mean_st(dir_ECCO)
    # Mask the missing data (zero values)
    SALT_extended[SALT_extended == 0.0] = np.nan
    THETA_extended[THETA_extended == 0.0] = np.nan

#change the order from (k_dep,i_lat,j_lon) to (j_lon,i_lat,k_dep),largely because
#gcc uses (lon,lat)
if not meanST_exist:
    SALT_tmp=np.transpose(SALT_mean[0,:,:,:],(2,1,0))
    THETA_tmp=np.transpose(THETA_mean[0,:,:,:],(2,1,0))
else:
    SALT_tmp=np.transpose(SALT_extended[:,:,:],(2,1,0))
    THETA_tmp=np.transpose(THETA_extended[:,:,:],(2,1,0))

#fill edges to extend the covered region.
print("Fill edges...")
SALT_lonlatdep=filledges(SALT_tmp)
THETA_lonlatdep=filledges(THETA_tmp)

# Compute great-circle path
print("Calculate geodetic path...")
d0, grid_path, lat_path, lon_path, x_path, p1_extend, p2_extend = compute_great_circle_path(psrc, prec, extend_m=extend_m, path_dx_m=path_dx_m)
path_xnpts=len(x_path)

print("Cut out a region from the globe to speed up computation...")
#cut out the relevant region from the globe, that helps speed up the calculation.
dlat_ECCO=lat_ECCO[1]-lat_ECCO[0]
dlon_ECCO=lon_ECCO[1]-lon_ECCO[0]
ilat_cut = np.where((lat_ECCO<=max(p1_extend[1],p2_extend[1])+2*dlat_ECCO) & (lat_ECCO>=min(p1_extend[1],p2_extend[1])-2*dlat_ECCO))[0]
ilon_cut = np.where((lon_ECCO<=max(p1_extend[0],p2_extend[0])+2*dlon_ECCO) & (lon_ECCO>=min(p1_extend[0],p2_extend[0])-2*dlon_ECCO))[0]

S_cut_ECCO = SALT_lonlatdep[ilon_cut[0]:ilon_cut[-1]+1,ilat_cut[0]:ilat_cut[-1]+1,:]
T_cut_ECCO = THETA_lonlatdep[ilon_cut[0]:ilon_cut[-1]+1,ilat_cut[0]:ilat_cut[-1]+1,:]
lat_cut_ECCO=lat_ECCO[ilat_cut[0]:ilat_cut[-1]+1]
lon_cut_ECCO=lon_ECCO[ilon_cut[0]:ilon_cut[-1]+1]

#use the gsw toolbox to compute sound speed and density
print("Call gsw to compute sound speed and density...")
c_cut_ECCO = c(T_cut_ECCO, S_cut_ECCO,lon_cut_ECCO.reshape(-1,1,1),lat_cut_ECCO.reshape(1,-1,1),-dep_ECCO)
rho_cut_ECCO = rho(T_cut_ECCO, S_cut_ECCO,lon_cut_ECCO.reshape(-1,1,1),lat_cut_ECCO.reshape(1,-1,1),-dep_ECCO)

########################find c,rho,theta, and salinity along the selected path###########
# --- Depth Interpolation ---
print("Do interpolation to get the profile of properties along the path...")
depth_path, c_path, rho_path, theta_path, salt_path = interpolate_properties_along_path(
    lon_cut_ECCO, lat_cut_ECCO, c_cut_ECCO, rho_cut_ECCO, T_cut_ECCO, S_cut_ECCO, 
    grid_path, dep_ECCO, fillbtm2, path_depth0_m, path_depth1_m, path_ddepth_m)
ndepth_path = depth_path.shape[0]
########################Done - find c,rho,theta, and salinity along the selected path###########


print("Save results...")
# --- Save Results ---
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, "path_results.npz")

#np.savez(output_file, 
#         sound_speed=c_pathx, 
#         density=rho_pathx, 
#         temperature=theta_pathx, 
#         salinity=salt_pathx, 
#         lat_path=lat_path, 
#         lon_path=lon_path, 
#         depths=dep_ECCO)
#print(f"Results saved to {output_file}")

speed_average, density_average, salt_average, theta_average = save_results(
    output_dir,
    tomo_file_name, 
    grid_file_name,
    ave1D_file_name, 
    grid_path,
    x_path, 
    depth_path, 
    c_path, 
    rho_path, 
    salt_path, 
    theta_path, 
    path_dx_m, 
    path_ddepth_m)

print("Sound speed, density, temperature, and salinity interpolation complete.")
