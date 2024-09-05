# load_pdbs_by_number.tcl

# List of folder numbers
set folder_numbers {84 159 190 195 242 252 288 316 339 371 450 477 518 537 544 632 638 655 879 888 893 905 908 930 953}

# Loop over each number and generate the folder name
foreach num $folder_numbers {
    # Generate the folder name based on the number
    set pdb_dir "./folder_trajectories_${num}_vmd"
    puts "Processing directory: $pdb_dir"

    # Change to the current directory
    cd $pdb_dir

    # Get the list of PDB files in the current directory
    set file_list [glob *.pdb]

    # Load each PDB file into VMD
    foreach file $file_list {
        puts "Loading file: $file"
        mol new $file type pdb
        set mol_id [molinfo top get id]

        # Set CPK representation for the molecule
        mol representation CPK
        mol color Name
        mol selection all
        mol addrep $mol_id
    }

    # Return to the parent directory
    cd ..
}

set clus_dir "./folder_MSM"
cd $clus_dir
set file_list [glob *.pdb]
foreach file $file_list {
    puts "Loading file: $file"
    mol new $file type pdb
    set mol_id [molinfo top get id]

    # Set CPK representation for the molecule
    mol representation CPK 3.0
    mol color ColorID 3
    mol selection all
    mol addrep $mol_id
}

set p1_dir "../"

# Change to the specified directory
cd $p1_dir

# Load the molecule
mol new p1_noh.pdb

# First representation: standard lines
mol representation Lines
mol color Name
set mol_id [molinfo top get id]
mol addrep $mol_id

# Second representation: newcartoon
mol representation NewCartoon
mol color Name
mol addrep $mol_id

# Third representation: CPK for specific residues
mol representation CPK 1.0 0.3 12.0 12.0
mol selection "resid 171 172 191 196"
mol color ColorID 4  ;# Yellow
mol addrep $mol_id
