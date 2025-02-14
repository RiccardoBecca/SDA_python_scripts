#!/bin/bash

# Loop over all sdamm_crowd_N.in files
for input_file in sdamm_crowd_*.in; do
    # Extract the number N from the filename
    N=$(echo "$input_file" | grep -oP '\d+')
    
    # Define the corresponding trajectory file
    trajectory_file="trajectories_$N"
    
    # Check if the corresponding trajectory file exists
    if [[ -f "$trajectory_file" ]]; then
        echo "Processing: $input_file and $trajectory_file"
        
        # Run the Python script
        python Monitor_contact.py \
            --sda_input_file "$input_file" \
            --trajectory "$trajectory_file" \
            --pdb_molecule p2_noh.pdb \
            --pdb_crowders crowder_noh.pdb \
            --mol_number 1 \
            --dist_coms 30.0 \
            --contact_dist 4.5 \
            --output_folder crowder_contacts
    else
        echo "Warning: Trajectory file $trajectory_file not found for input $input_file. Skipping."
    fi
done
