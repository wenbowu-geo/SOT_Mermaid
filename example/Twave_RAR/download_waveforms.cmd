awk 'BEGIN{start_write=0;} { if(start_write==1) {split($3,evt_time,",");  split($0,tmp,","); lat=tmp[6];lon=tmp[7]; evdp=tmp[8]; split($0,mag_tmp,"mb"); split(mag_tmp[2],mag,","); if(mag[2]=="") mag[2]=0.0; print evt_time[2],evt_time[3],lat,lon,evdp,mag[2];} if($1=="EVENTID,TYPE,AUTHOR") start_write=1;   }' ISC_catalog.txt >ISC_catalog_reformat.txt

python download_data_given_dist_stations.py
ls waveforms/*/*mseed >mseed_list.txt
python mseedTOsac.py mseed_list.txt

#add earthquake source parameters to sac header.
awk '{split($1,part1,"-"); split($2,part2,":");  folder="waveforms/"part1[1]""part1[2]""part1[3]""part2[1]""part2[2]""part2[3]; print "sh $download_data_obspy//add_evt_info/add_oneevt.cmd", $1,$2,$4,$3,$5,$6,folder;}' ISC_catalog_reformat.txt |sh

#add station lon/lat to sac header.
ls -d waveforms/*  >dir_list
awk '{print "ls", $1"/*?H?",">sac_trace_list"; print "python $download_data_obspy/add_stla_stlo/add_stla_stlo.py"}' dir_list |sh
#set lcalda with true and sac would automatically compute gcarc and az using (evla,evlo,stla,stlo)
awk '{print "ls", $1"/*?H?";}' dir_list |sh | awk '{print "r ",$1; print "ch lcalda true"; print "w over";}' |sac

#mark T1 as T-wave arrival time 
saclst dist f waveforms/*/*HZ >distance.txt
awk '{print "r",$1; print "ch t1",$2/1.5; print "w over";}' distance.txt |sac
