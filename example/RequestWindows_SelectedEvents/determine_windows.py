from requested_time_windows import requested_time_windows
import os

#sound speed to predict T-arrival
speed=1.5

#The below time window is selected for requesting.
DTbefore=70.0 #second before T-wave arrival
DTafter=170.0 #second after T-wave arrival

#Floats in the below file are selected for requesting. 
floats_list_file='data/floats.txt'

#earthquakes in the below catalog are selected for requesting.
eq_catalog_file='data/ISC_catalog.txt'


#####Ensure the Mermaid trajectory files, bathymetry .grd file (used for plotting maps) and 
#####Mermaid acquization file are ready under the the below path.
path_global_data='../../data/'


#Path of results
path_results='results/'
path_figure='results/figures/'
os.makedirs(path_results, exist_ok=True)
os.makedirs(path_figure, exist_ok=True)

#plot maps or not.
plot_map=False

#####Find the time windows.
time_windows = requested_time_windows(speed,DTbefore,DTafter,floats_list_file,
                  eq_catalog_file,path_global_data,plot_map,path_figure)

# Save the requested time windows
file_request_window = 'results/request_window.txt'
with open(file_request_window, 'w') as file_request_window_handle:
    file_request_window_handle.write('Float TimeBegin TimeStop eq_id\n')
    for float_key, timestamp_list in time_windows.items():
        for tstart, tend, ieq in timestamp_list:
            file_request_window_handle.write('{} {} {} {}\n'.format(
                float_key,
                tstart.strftime('%Y-%m-%dT%H:%M:%S'),
                tend.strftime('%Y-%m-%dT%H:%M:%S'),ieq
            ))
