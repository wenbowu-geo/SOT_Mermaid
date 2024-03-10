ls P*_geo_DET_REQ.csv > filenames.txt
matlab -nosplash -nodesktop -r "run('extract_trajectory.m'); exit;"

