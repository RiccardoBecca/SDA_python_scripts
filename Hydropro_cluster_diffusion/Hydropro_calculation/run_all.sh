#!/bin/bash

# Define the main directory where Hydropro_calculation is located
#MAIN_DIR="$(dirname "$(realpath "$0")")"

# Loop over directories matching the pattern *bounded*
for folder in "."/*_bounded; do
    # Ensure it's a directory
    if [[ -d "$folder" ]]; then
        echo "Processing folder: $folder"
        
        # Loop over all PDB files in the current bounded folder
        for pdb_file in "$folder"/*.pdb; do
            # Ensure there are PDB files
            if [[ -f "$pdb_file" ]]; then
                echo "Processing file: $pdb_file"

                # Copy the PDB file to Hydropro_calculation and rename it to tmp.pdb
                cp "$pdb_file" "./Hydropro_calculation/tmp.pdb"

                # Navigate to Hydropro_calculation and run the script
                (cd "./Hydropro_calculation" && ./run_hydropro.sh tmp.pdb)
		python ./Hydropro_calculation/extract_diff_coeff.py --output_file tmp-sda_values.dat --folder_name "./diffusion_coefficients/$folder" --file_name "./diffusion_coefficients/$pdb_file"
            fi
        done
    fi
done

cd "./Hydropro_calculation"
rm hydropro.dat hydropro-fit.txt hydropro-sum.txt tmp*

echo "All done!"

