#!/bin/bash

# Input files and directories
file_sta="data/loc.txt"
dir_ECCO="../../data/ECCO/"
file_sta="data/loc_seafloor_cutoff.txt"
output_dir="results/"

# Ensure output directory exists and clean previous results
mkdir -p "${output_dir}"
rm -f "${output_dir}/average1D_VpRhoTS.txt"	
rm -f "${output_dir}/mode_1_*Hz.txt"	
rm -f "${output_dir}/tomo_acoustic.txt"

# Define frequency parameters
nfreqs=4
freq_min=2.5
freq_max=10.0

# Define sediment structure (depth, velocity, density)
thickness_sedi=100.0  # m
vp_sedi=5800.0        # m/s
vs_sedi=3200.0        # m/s
rho_sedi=2.6          # g/cm^3

# Define lower half-space structure
vp_lowerhalf=5800.1   # m/s 
vs_lowerhalf=3200.1   # m/s
rho_lowerhalf=2.61     # g/cm^3

# Distance (in meters) to extend the path range north and south from each station
path_extend_distance=50000.0 #unit of meter

# Seafloor cutoff depth setting
seafloor_cutoff_etopo1=False

# If using user-specified cutoff depth, ensure the cutoff file is generated
if [ "$seafloor_cutoff_etopo1" = False ]; then
    echo "Using ETOPO1 for seafloor cutoff depth."
    if [ ! -f "$file_sta" ]; then
        echo "Error: $file_sta not found! Run cutoff depth extraction first."
        exit 1
    fi
else
    echo "Generating seafloor cutoff depth using ETOPO1..."
    echo "Reading data/loc_Mermaid.txt to generate data/loc_seafloor_cutoff.txt"
    python ../../src/cutoff_depth_etopo1.py --dir_etopo1 "../../data/Bahymetry/" --extend_m ${path_extend_distance}
fi
echo

# Initialize pair counter
ipair=0

# Loop over each line in loc.txt, skipping the header.
tail -n +2 ${file_sta} | while IFS=" " read -r sta lat_rec lon_rec dep_rec evl_rec; do
    # Increment the source coordinates for each pair (example adjustment).
    lon_src=$(echo "$lon_rec + 0.0" | bc)
    lat_src=$(echo "$lat_rec + 0.05" | bc)
    echo ${lat_rec} ${evl_rec}

    # Increment the pair counter
    ipair=$((ipair + 1))
    echo "Running pair $ipair: lon_src=$lon_src lat_src=$lat_src lon_rec=$lon_rec lat_rec=$lat_rec"

    # Set the cutoff depth to the seafloor depth
    dep_cutoff=$(echo "-1 * $evl_rec" | bc)
    #Run Python script to generate sound speed profile
    python ../../src/SoundSpeed_ECCO.py --psrc "$lon_src" "$lat_src" --prec "$lon_rec" "$lat_rec" \
               --extend_m ${path_extend_distance} --path_dx_m 5000.0 \
               --path_depth0_m 0.0 --path_depth1_m ${dep_cutoff} --path_ddepth_m 10.0 \
               --meanST_exist True --year0_mean 1992 --year1_mean 1992 \
               --dir_ECCO "${dir_ECCO}" --output_dir "${output_dir}"

    # Change to results/ directory
    pushd "${output_dir}" > /dev/null

    # Run MATLAB script with properly formatted arguments
    matlab -nodisplay -r "addpath('/Users/wenbowu/work/T_wave/Mermaid/SOT_Mermaid/src/'); \
         ORCA_DeepOcean_Modes(${nfreqs}, ${freq_min}, ${freq_max}, ${thickness_sedi}, \
         ${vp_sedi}, ${vs_sedi}, ${rho_sedi}, ${vp_lowerhalf}, ${vs_lowerhalf}, ${rho_lowerhalf}); exit"


    # Create directory for results
    dir_name="${sta}"
    mkdir -p "${dir_name}"

    # Move generated files into the station directory
    mv mode_*txt "${dir_name}/"
    mv tomo_acoustic.txt "${dir_name}/"
    mv average1D_VpRhoTS.txt "${dir_name}/"
    mv property_path.npz "${dir_name}/"

    # Return to previous directory
    popd > /dev/null
done

