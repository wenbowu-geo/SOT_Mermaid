#!/bin/bash

# Input files and directories
file_sta="data/loc.txt"
dir_ECCO="../../data/ECCO/"
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
vp_sedi=5600.0        # m/s
vs_sedi=3200.0        # m/s
rho_sedi=2.6          # g/cm^3

# Define lower half-space structure
vp_lowerhalf=5600.1   # m/s 
vs_lowerhalf=3200.1   # m/s
rho_lowerhalf=2.61     # g/cm^3

# Initialize pair counter
ipair=0

# Process each line in loc.txt, skipping the header
tail -n +2 "${file_sta}" | while IFS=' ' read -ra line; do
    # Ensure the line has enough columns
    if [ "${#line[@]}" -lt 4 ]; then
        echo "Skipping malformed line: ${line[*]}"
        continue
    fi

    # Read station information
    sta="${line[0]}"
    lat_rec="${line[1]}"
    lon_rec="${line[2]}"
    dep_rec="${line[3]}"

    # Adjust source location slightly (use awk for floating-point math)
    lon_src=$(awk -v lon="$lon_rec" 'BEGIN { printf "%.6f", lon + 0.0 }')
    lat_src=$(awk -v lat="$lat_rec" 'BEGIN { printf "%.6f", lat + 0.05 }')

    # Increment pair counter
    ipair=$((ipair + 1))
    echo "Running pair $ipair: lon_src=$lon_src lat_src=$lat_src lon_rec=$lon_rec lat_rec=$lat_rec"

    #Run Python script to generate sound speed profile
    python ../../src/SoundSpeed_ECCO.py --psrc "$lon_src" "$lat_src" --prec "$lon_rec" "$lat_rec" \
               --extend_m 50000.0 --path_dx_m 5000.0 \
               --path_depth0_m 0.0 --path_depth1_m 8000.0 --path_ddepth_m 10.0 \
               --meanST_exist True --year0_mean 1992 --year1_mean 1992 \
               --dir_ECCO "${dir_ECCO}" --output_dir "${output_dir}"

    # Change to results/ directory
    pushd "${output_dir}" > /dev/null

    # Run MATLAB script with properly formatted arguments
    #matlab -nodisplay -r "addpath('/Users/wenbowu/work/T_wave/Mermaid/SOT_Mermaid/src/'); \
    #    DeepOcean_Modes(${nfreqs}, ${freq_min}, ${freq_max}, ${thickness_sedi}, \
    #    ${vp_sedi}, ${vs_sedi}, ${rho_sedi}, ${vp_lowerhalf}, ${vs_lowerhalf}, ${rho_lowerhalf}); exit;"
    matlab -nodisplay -r "addpath('/Users/wenbowu/work/T_wave/Mermaid/SOT_Mermaid/src/'); \
         DeepOcean_Modes(${nfreqs}, ${freq_min}, ${freq_max}, ${thickness_sedi}, \
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

