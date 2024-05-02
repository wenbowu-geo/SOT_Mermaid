import pandas as pd
from datetime import timedelta
import numpy as np
from geopy.distance import geodesic
from read_etopo1 import read_etopo1_netcdf
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import copy
import matplotlib.colors as colors

def requested_time_windows(speed,DTbefore,DTafter,floats_list_file,eq_catalog_file,path_global_data,plot_map,path_figure):
    """
    Determine requested time windows of Mermaid floats for T-waves.

    Parameters:
    :param float speed: Sound speed, in units of km/s.
    :param float DTbefore: Time before the predicted T-wave arrival for requesting data, in units of seconds.
    :param float DTafter: Time after the predicted T-wave arrival for requesting data, in units of seconds.
    :param str floats_list_file: File name in which the selected floats' names are listed.
    :param str eq_catalog_file: File name of the earthquake catalog. The catalog has the format of ISC report.
                                The first 25 lines are header and skipped in reading.
    :param str path_global_data: The path name of the globally used database. This function requires bathymetry data,
                                 Mermaid trajectories, and acquisition log files, all of which are saved under the path
                                 path_global_data.
    :param bool plot_map: If True, plot raypath maps.
    :param str path_figure: The path name to save the raypath maps.

    Returns:
    :return dict request_windows: Dictionary to save the requested time windows. Its structure is like
                                   request_windows['P0010'][ieq][0:2], where the key 'P0010' is the float name,
                                   the index 'ieq' is the ith earthquake, and the last dimension corresponds to
                                   start datetime and end datetime of the requested window and the equauake's index 
                                   in the catalog.

    Usage:
    requested_time_windows = requested_time_windows(1.5, 70.0, 170.0, 'data/floats.txt', 'data/ISC_catalog.txt',
                                                     '../../data/', True, 'results/figures/')
    """


    ####################################Read target floats#######################
    floats = np.atleast_1d(np.genfromtxt(floats_list_file, dtype=str))
    
    ####################################Read earthquake catalogs#################
    with open(eq_catalog_file, 'r') as file:
        column_names = file.readlines()[25].strip().split(',')
    column_names = [name.strip() for name in column_names]
    column_names[-2]="TYPE_MAG"
    column_names[-3]="AUTHOR_MAG"
    # Read the ISC catalog into a DataFrame, skipping the first 25 rows
    eq = pd.read_csv(eq_catalog_file, skiprows=26, header=None, names=column_names,usecols=range(len(column_names)))
    # Convert 'date' and 'time' columns in eq to datetime format
    eq['datetime'] = pd.to_datetime(eq['DATE'] + ' ' + eq['TIME'], format='%Y-%m-%d %H:%M:%S.%f')
    eq.drop(['DATE', 'TIME'], axis=1, inplace=True)
    
    ####################################Read Mermaid trajectories################
    traj = {}
    for key in floats:
        file = path_global_data + '/Mermaid_coordinates/' + key + '_traj.txt'
        traj[key] = pd.read_csv(file, delim_whitespace=True, skiprows=1, names=['date', 'lat', 'lon', 'dep'])
        traj[key]['datetime'] = pd.to_datetime(traj[key]['date'], format='%Y-%m-%dT%H:%M:%S.%fZ')
    
    ###################################Determine the requested time wondows######
    # Create an empty nested dictionary to store the data availability.
    data_missing = {}
    # Create an empty nested dictionary to store the requested windows.
    request_windows={}
    float_loc={}
    for float_key in floats:
        print('float:',float_key)
        # Do interpolation to find float's lat/lon at the earthquake's origin time
        interp_lat = np.interp(eq['datetime'],traj[float_key]['datetime'],traj[float_key]['lat'])
        interp_lon = np.interp(eq['datetime'],traj[float_key]['datetime'],traj[float_key]['lon'])
        float_loc[float_key] = np.column_stack((interp_lon, interp_lat))
        interp_coordinates = list(zip(interp_lat, interp_lon))
        # Read this float's acquisition log file
        file_path = path_global_data + '/Mermaid_act_log/' + float_key + '_acq_log.txt' 
        with open(file_path, 'r') as file:
            lines = file.readlines()
        acq = pd.DataFrame([(line.split() + [None, None])[:2] for line in lines], columns=['Timestamp', 'Event'])
        # Convert the 'Timestamp' column to datetime
        acq['Timestamp'] = pd.to_datetime(acq['Timestamp'].str[:19], format='%Y-%m-%dT%H:%M:%S')
        def adjust_timestamp(row):
            if row['Event'] == 'started':
                return row['Timestamp'] + timedelta(minutes=2)
            else:
                return row['Timestamp'] 
        # Apply the 2-min delay adjustment to 'started'
        acq['Timestamp'] = acq.apply(adjust_timestamp, axis=1)
        # Create a list to store available data sections
        data_sections = []
        # Initialize variables to track the start and end of a section
        section_start = None; 
        # Track the current status of recording. 
        # Occasionally, two consequtive 'started' appear, and the first will be ignored. 
        start_recording=False
        # Iterate through the DataFrame to create data sections
        for index, row in acq.iterrows():
            timestamp, event = row['Timestamp'], row['Event']
            if event == 'started':
                section_start = timestamp
                start_recording=True
            elif start_recording:
                section_end = timestamp
                data_sections.append([section_start, section_end])
        # Create a list to store data availability and requested windows for this float
        data_missing_this_float = []; windows_this_float=[]
        for ieq in range(len(eq['datetime'])):
            distance=geodesic((eq['LAT'][ieq],eq['LON'][ieq]),(interp_lat[ieq],interp_lon[ieq])).km
            Tarrival=eq['datetime'][ieq]+pd.to_timedelta(distance/speed, unit='s')
            tstart=Tarrival-pd.to_timedelta(DTbefore, unit='s')
            tend=Tarrival+pd.to_timedelta(DTafter, unit='s')
            # Check if the selected window fully falls within any data section
            window_fully_in_data_section = any((tstart >= start) and (tend <= end) for start, end in data_sections)
            # Update the data availability list for this float
            data_missing_this_float.append(not window_fully_in_data_section)
            print('    earthquake:', eq['datetime'][ieq])
            print('    T-wave arrival:',eq['datetime'][ieq]+pd.to_timedelta(distance/speed, unit='s'))
            print('    distance={:.2f} km, T-wave arrival={:.2f} s.'.format(distance, distance/speed))
            print('    requested window: {} to {}, data missing: {}.'.format(tstart.strftime('%Y-%m-%dT%H:%M:%S'),
                   tend.strftime('%Y-%m-%dT%H:%M:%S'),not window_fully_in_data_section))
            print()
            # Save the requsted time windows
            if window_fully_in_data_section:
                 windows_this_float.append([tstart,tend,ieq])
        # Add this float's data availability and requested time window lists to the corresponding dictionaries.
        data_missing[float_key] = data_missing_this_float
        request_windows[float_key] = windows_this_float
        print()
    
    ###########################################plot maps##############################
    if plot_map:
        def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
            new_cmap = colors.LinearSegmentedColormap.from_list(
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                cmap(np.linspace(minval, maxval, n)))
            return new_cmap
        
        fontsize=7
        plt.rcParams.update({'font.size': fontsize})
        
        etopo_file= path_global_data + "/Bahymetry/etopo5.grd"
        #-180 to 180 deg if lon0>180, or 0-360 deg if lon0<=180.0
        lon0=170
        lon_etopo1,lat_etopo1,ele_etopo1=read_etopo1_netcdf(etopo_file,170)
        dlon_etopo1=lon_etopo1[1]-lon_etopo1[0]
        dlat_etopo1=lat_etopo1[1]-lat_etopo1[0]

        #use -180 to 180 deg for lon.
        lon_etopo1_shifted = lon_etopo1 - 180.0
        lon_etopo1_shifted[lon_etopo1_shifted < -180.0] += 360.0
        etopo_lon_grid, etopo_lat_grid = np.meshgrid(lon_etopo1_shifted, lat_etopo1)
        
        #cmap = plt.get_cmap('Greys')
        #new_cmap = truncate_colormap(cmap, 0., 0.5)
        orig_map=plt.cm.get_cmap('Blues')
        reversed_map=orig_map.reversed()
        
        for ieq in range(len(eq['datetime'])):
            print(ieq)
            projection = ccrs.PlateCarree(central_longitude=180)
            fig, ax0 = plt.subplots(subplot_kw={'projection': projection})
            im = ax0.pcolormesh(etopo_lon_grid,etopo_lat_grid, ele_etopo1, vmin=-8000.0, vmax=1000., cmap=reversed_map, shading='auto')
            #ax0.coastlines(resolution='10m', color='black', linewidth=1)
            #ax0.add_feature(cartopy.feature.BORDERS, linestyle=':', linewidth=1)
            x_evt, y_evt = ax0.projection.transform_point(eq['LON'][ieq], eq['LAT'][ieq], ccrs.PlateCarree())
            ax0.scatter(x_evt,y_evt, s=100, marker="o", color='black')
            # Draw meridians and parallels
            ax0.set_xticks(np.arange(-180, 181, 10), crs=ccrs.PlateCarree())
            ax0.set_yticks(np.arange(-90, 91, 5), crs=ccrs.PlateCarree())
            ax0.xaxis.set_major_formatter(LongitudeFormatter())
            ax0.yaxis.set_major_formatter(LatitudeFormatter())
            for float_key in floats:
                x, y = ax0.projection.transform_point(float_loc[float_key][ieq][0], float_loc[float_key][ieq][1], ccrs.PlateCarree())
                ax0.scatter(x, y, s=100, color='black', marker="^")
                if data_missing[float_key][ieq]:
                    ax0.plot([x, x_evt], [y, y_evt], color='tab:red', linewidth=1.0, label='data gap')
                else:
                    ax0.plot([x, x_evt], [y, y_evt], color='black', linewidth=1.0, label='available')
                ax0.text(x + 1, y, float_key, color='black', fontsize=fontsize)
            float_lons = [float_loc[key][ieq][0] for key in float_loc]
            float_lats = [float_loc[key][ieq][1] for key in float_loc]
            lons_rec_src=float_lons+[eq['LON'][ieq]]
            #adjust to 0-360 deg range.
            lons_rec_src = [(lon + 360) % 360 for lon in lons_rec_src]
            lats_rec_src=float_lats+[eq['LAT'][ieq]]
            ax0.set_extent([np.min(lons_rec_src)-5, np.max(lons_rec_src)+5, np.min(lats_rec_src)-5,np.max(lats_rec_src)+5], ccrs.PlateCarree())
            ax0.set_title(f"{eq['datetime'][ieq]} , Magnitude: {eq['MAG'][ieq]}, red raypath for missing data", color='black', fontsize=fontsize)
            plt.savefig(path_figure + f"{eq['datetime'][ieq].strftime('%Y-%m-%dT%H:%M:%S.%f')}.png", dpi=400, format='png')
            #plt.show()
            plt.close(fig)
    return request_windows
