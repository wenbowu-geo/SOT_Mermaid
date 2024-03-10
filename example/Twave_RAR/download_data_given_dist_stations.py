#!/usr/bin/env python
import obspy
from obspy.clients.fdsn.client import Client
from obspy.clients.fdsn.mass_downloader import CircularDomain, \
    Restrictions, MassDownloader

import os
import warnings
import datetime

def data_download(stations, starttime, endtime, event_name):

    print("\n========================================")
    print("event:", event_name)
    print("time:", starttime, endtime)
    waveforms_folder = "waveforms/" + event_name
    stationxml_folder = "stations/" + event_name
    c = Client("IRIS")
    #clients could be one of IRIS, GFZ, ORFEUS, ETH, GEONET etc.

    if not os.path.exists(waveforms_folder):
        os.makedirs(waveforms_folder)

    if not os.path.exists(stationxml_folder):
        os.makedirs(stationxml_folder)

    # First download waveforms.
    for network, station in stations:
        filename = os.path.join(waveforms_folder,
                            "%s.%s.mseed" % (network, station))
        if os.path.exists(filename):
            continue

        try:
            c.get_waveforms(network=network, station=station, location="*",
                            channel="EDH", starttime=starttime, endtime=endtime,
                            filename=filename)
        except Exception as e:
            print("Failed to download %s.%s due to %s" %
                (network, station, str(e)))
            continue

        print("Successfully downloaded %s." % filename)

        stationxml_filename = os.path.join(stationxml_folder,
                                       "%s.%s.xml" % (network, station))

        if os.path.exists(stationxml_filename):
            continue

        try:
            c.get_stations(network=network, station=station, location="*",
                            channel="?H?,EDH", starttime=starttime, endtime=endtime,
                            filename=stationxml_filename, level="response")
        except Exception as e:
            print("Failed to download %s.%s StationXML due to %s" % (
                network, station, str(e)))
            continue

        print("Successfully downloaded %s." % stationxml_filename)


def download_given_dist_sta(evt_lat,evt_lon,start_time,end_time,event_name,\
                             given_min_dist,given_max_dist,given_network,given_station):

# Circular domain around the epicenter. This will download all data between
# 70 and 90 degrees distance from the epicenter. This module also offers
# rectangular and global domains. More complex domains can be defined by
# inheriting from the Domain class.
    domain = CircularDomain(latitude=evt_lat, longitude=evt_lon, 
                        minradius=given_min_dist, maxradius=given_max_dist)

    restrictions = Restrictions(
    # Get data from 5 minutes before the event to one hour after the
    # event. This defines the temporal bounds of the waveform data.
    #starttime=origin_time - 5 * 60,
    #endtime=origin_time,
    starttime=start_time,
    endtime=end_time+7000,
    #network="IM,GE,IU,AU,II,PS,MY",
    #station='H09W1,H09N1,H08S1,H08S2,H08S3,H08N1,H08N2,H08N3,H01W1,H01W2,H01W3,COCO,XMIS,DGAR,PALK,HMDM,KAAM,HALK,MALK,GSI,MNC,PSI,KUM',
    network=given_network,
    station=given_station,
    # You might not want to deal with gaps in the data. If this setting is
    # True, any trace with a gap/overlap will be discarded.
     reject_channels_with_gaps=True,
    # And you might only want waveforms that have data for at least 95 % of
    # the requested time span. Any trace that is shorter than 95 % of the
    # desired total duration will be discarded.
     minimum_length=0.95,
    # No two stations should be closer than 10 km to each other. This is
    # useful to for example filter out stations that are part of different
    # networks but at the same physical station. Settings this option to
    # zero or None will disable that filtering.
     minimum_interstation_distance_in_m=1E3,
    # Only HH or BH channels. If a station has HH channels, those will be
    # downloaded, otherwise the BH. Nothing will be downloaded if it has
    # neither. You can add more/less patterns if you like.
     #channel_priorities=["BH?","HH?"],
     channel_priorities=["BH?"],
    # Location codes are arbitrary and there is no rule as to which
    # location is best. Same logic as for the previous setting.
     location_priorities=["", "00", "10"])

    mseed_storage = "waveforms/"+event_name
# No specified providers will result in all known ones being queried.
    mdl = MassDownloader(providers=["IRIS"])
    #mdl = MassDownloader(providers=["ORFEUS"])
# The data will be downloaded to the ``./waveforms/`` and ``./stations/``
# folders with automatically chosen file names.
    mdl.download(domain, restrictions, mseed_storage="waveforms/"+event_name,
             stationxml_storage="stations/"+event_name)


def read_station_file(station_filename):
    stations = []
    with open(station_filename, "rt") as fh:
        for line in fh:
            line = line.split()
            stations.append((line[1], line[0]))
    return stations


def read_cmt_file(filename):
    """
    Initialize a source object from a CMTSOLUTION file.
    :param filename: path to the CMTSOLUTION file
    """
    with open(filename, "rt") as f:
        line = f.readline()
        origin_time = line[4:].strip().split()[:6]
        values = list(map(int, origin_time[:-1])) + [ float(origin_time[-1]) ]
        try:
            origin_time = obspy.UTCDateTime(*values)
            print(origin_time)
        except (TypeError, ValueError):
            warnings.warn("Could not determine origin time from line: %s" % line)
            origin_time = obspy.UTCDateTime(0)
        event_name = f.readline().split()[-1]
        # print type(data), data
        #event_name = float(f.readline().strip().split()[-1])
        time_shift = float(f.readline().strip().split()[-1])
        half_duration = float(f.readline().strip().split()[-1])
        latitude = float(f.readline().strip().split()[-1])
        longitude = float(f.readline().strip().split()[-1])
        depth_in_m = float(f.readline().strip().split()[-1]) * 1e3
        m_rr = float(f.readline().strip().split()[-1]) / 1e7
        m_tt = float(f.readline().strip().split()[-1]) / 1e7
        m_pp = float(f.readline().strip().split()[-1]) / 1e7
        m_rt = float(f.readline().strip().split()[-1]) / 1e7
        m_rp = float(f.readline().strip().split()[-1]) / 1e7
        m_tp = float(f.readline().strip().split()[-1]) / 1e7
    return (latitude, longitude, depth_in_m, m_rr, m_tt, m_pp, m_rt,
            m_rp, m_tp, time_shift, origin_time, event_name, half_duration)


def read_cmt_list(cmt_list_file):
    cmt_list = []
    with open(cmt_list_file) as f:
        content = f.read().splitlines()
        #print "content:", content
        #cmt_list.append(content)
    return content

#WENBO
def read_evtcata_list(evtcat_file):
    evt_list = []
    with open(evtcat_file) as f:
        content = f.read().splitlines()
        #print "content:", content
        #cmt_list.append(content)
    return content

#station_filename = "./STATIONS.txt"
#station_list = read_station_file(station_filename)
#print "station_list:", station_list

#cmt_list_file = "./cmt_list"
#cmt_list = read_cmt_list(cmt_list_file)
#print "cmt_info:", cmt_list

#for cmt in cmt_list:
#    cmt_file = './cmt/' + cmt
#    cmt_info = read_cmt_file(cmt_file)
#    event_name = cmt_info[11]
#    calibrate_time = cmt_info[10] + datetime.timedelta(seconds=int(cmt_info[9]))
#    starttime = calibrate_time + datetime.timedelta(seconds=-20)
#    endtime = calibrate_time + datetime.timedelta(seconds=120)
#    data_download(station_list, starttime, endtime, event_name)

evt_list_file = "./ISC_catalog_reformat.txt"
evt_list = read_evtcata_list(evt_list_file)
print("evt_info:",evt_list_file)

for evt in evt_list:
     info=evt.split()
     values = info[0]+'T'+info[1]
     try:
        origin_time = obspy.UTCDateTime(values)
        print(origin_time)
     except (TypeError, ValueError):
        warnings.warn("Could not determine origin time from line: %s" % line)
        origin_time = obspy.UTCDateTime(0)
     yeardate=info[0].split("-")
     time=info[1].split(":")
     event_name=yeardate[0]+yeardate[1]+yeardate[2]+time[0]+time[1]+time[2]
#     print event_name
     starttime = origin_time + datetime.timedelta(seconds=-5)
     endtime = origin_time + datetime.timedelta(seconds=20)
     evt_lat=info[2]
     evt_lon=info[3]
#     data_download(station_list, starttime, endtime, event_name)

     #download the data of the below stations. Most of them, except PSI and KUM,
     # have the potential record T-waves. PSI and KUM are used to search repeaters.
     #given_network="X9,7D,II,IU,TA,US,CI,N4,PN,LD,ZC"
     given_network="IU"
     #given_station='OBS56,OBS57'
     given_station='RAR'
     given_min_dist=0.0
     given_max_dist=90.0
     download_given_dist_sta(evt_lat,evt_lon,starttime,endtime,event_name,\
                             given_min_dist,given_max_dist,given_network,given_station)

     #download the data of all the stations within a distance of 10 deg.
     #Some of these stations might be already downloaded previously and that is fine.
     given_network="*"
     given_station='*'
     given_min_dist=0.0
     given_max_dist=0.3
     #download_given_dist_sta(evt_lat,evt_lon,starttime,endtime,event_name,\
     #                        given_min_dist,given_max_dist,given_network,given_station)



#print cmt_info
#print starttime, endtime
