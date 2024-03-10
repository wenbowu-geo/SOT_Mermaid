import obspy
from obspy import read, read_inventory
import os
import sys


def read_trace_list(trace_list_file):
    trace_list = []
    with open(trace_list_file) as f:
        content = f.read().splitlines()
    return content

#trace_list_file = "./mseed_trace_list"
trace_list_file=sys.argv[1]
trace_list = read_trace_list(trace_list_file)
for trace in trace_list:
     tracesplit_folder=trace.split("/")
     nlevels_path=len(tracesplit_folder)
     path=tracesplit_folder[0]+'/'
     for ilevel_path in range(1,nlevels_path-1):
         path=path+tracesplit_folder[ilevel_path]+'/'
     net_sta=tracesplit_folder[nlevels_path-1].split(".")
     try:
         #print 'work on event',tracesplit_folder[nlevels_path-2],',station',
         #       net_sta[1],' of network',net_sta[0]
         print('work on event',tracesplit_folder[nlevels_path-2],',station', \
               net_sta[1],',network',net_sta[0])

         st = read(trace)
         sta_response="stations/"+tracesplit_folder[1]+"/"+net_sta[0]+"."+net_sta[1]+".xml"
         for i_mseed in range(len(st)):
             tr_this_channel = st[i_mseed].copy()
             this_channel=tr_this_channel.stats.channel
             rename_trace=path+net_sta[0]+"."+net_sta[1]+".."+this_channel
             try:
                 inv = read_inventory(sta_response)
                 pre_filt = [0.001, 0.005, 10, 20]
                 #tr.remove_response(inventory=inv, pre_filt=pre_filt, output="VEL",
                 #      water_level=60, plot=False)  
                 tr_this_channel.remove_response(inventory=inv, output="VEL",
                      water_level=60, plot=False)  
                 tr_this_channel.write(rename_trace,format='SAC')
             except Exception:
                 print('failed, removing response:',rename_trace)
                 pass

     except Exception:
         print('failed, station',net_sta[1],' of network',net_sta[0])
         pass

