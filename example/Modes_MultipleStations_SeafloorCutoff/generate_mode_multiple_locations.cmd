#!/bin/bash
# Define an array of frequencies for mode computations.
freqs=(2.5 5.0 7.5 10.0)
file_sta="data/loc.txt"
dir_ECCO="../../data/ECCO/"

output_dir="results/"
mkdir -p ${output_dir}
rm -f ${output_dir}/average1D_VpRhoTS.txt	
rm -f ${output_dir}/mode_1_*Hz.txt	
rm -f ${output_dir}/tomo_acoustic.txt

#Lower half-space structure
Vp_lowerhalf=5800.0
Vs_lowerhalf=1000.0
Rho_lowerhalf=2.6

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
    # Run the Python script with the current location pair.
    python ../../src/SoundSpeed_ECCO.py --psrc "$lon_src" "$lat_src" --prec "$lon_rec" "$lat_rec" \
               --extend_m 50000.0 --path_dx_m 5000.0 \
               --path_depth0_m 0.0 --path_depth1_m ${dep_cutoff} --path_ddepth_m 10.0 \
               --meanST_exist True --year0_mean 1992 --year1_mean 1992 \
               --dir_ECCO ${dir_ECCO} \
               --output_dir ${output_dir}

    ######################## Prepare input file for running the MATLAB script KRAKEN #################
    # Write the first part of Twave_example.env
    echo "'Twave'" > ${output_dir}/Twave_example.env
    echo '2.0 #frequency' >> ${output_dir}/Twave_example.env
    echo '1' >> ${output_dir}/Twave_example.env
    echo "'NVF'" >> ${output_dir}/Twave_example.env

    # Get maximum depth from the processed input file
    max_dep=$(awk '{if(NR>2 && $3!="nan") printf("%.2f %.2f 0.0 1.0 /\n",$1, $2);}' ${output_dir}/average1D_VpRhoTS.txt | awk '{max_dep=$1;} END{print max_dep;}')
    echo "2000     0.000000 $max_dep" >> ${output_dir}/Twave_example.env

    # Append processed data from the input file
    awk '{if(NR>2 && $3!="nan") printf("%.2f %.2f 0.0 1.0 /\n",$1, $2);}' ${output_dir}/average1D_VpRhoTS.txt >> ${output_dir}/Twave_example.env

    # Append the final part of Twave_example.env
    echo "'A' 0.0" >> ${output_dir}/Twave_example.env
    echo "$max_dep $Vp_lowerhalf $Vs_lowerhalf $Rho_lowerhalf /" >> ${output_dir}/Twave_example.env
    echo '1300 2000' >> ${output_dir}/Twave_example.env
    echo '0' >> ${output_dir}/Twave_example.env
    echo '1' >> ${output_dir}/Twave_example.env
    echo '500 /' >> ${output_dir}/Twave_example.env
    echo '400' >> ${output_dir}/Twave_example.env
    echo '0 8000.0 /' >> ${output_dir}/Twave_example.env

    ######################## Run acoustic mode for each frequency ################
    #run KRAKEN_mode_func.m in the folder results/
    # Change to results/ directory to ensure KRAKEN_mode_func has access to required files
    pushd "$output_dir" > /dev/null
    for freq in "${freqs[@]}"; do
        # Run the MATLAB script with the current frequency
        matlab -nodisplay -r "addpath('/Users/wenbowu/work/T_wave/Mermaid/SOT_Mermaid/src/'); KRAKEN_mode_func($freq); exit;"
    done
    ######################## Create a directory and move files #################
    dir_name="$sta"
    mkdir -p "$dir_name"  # Create directory if it doesn't exist
    mv mode_*txt ${dir_name}/
    mv tomo_acoustic.txt ${dir_name}/
    mv average1D_VpRhoTS.txt ${dir_name}/
    mv property_path.npz ${dir_name}/
    popd > /dev/null
done
