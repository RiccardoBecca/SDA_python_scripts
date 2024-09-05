
set clus_dir "./folder_msm"
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
