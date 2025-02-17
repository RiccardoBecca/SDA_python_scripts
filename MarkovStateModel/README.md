# Markov State Model for SDA.



<img src="https://github.com/RiccardoBecca/SDA_python_scripts/blob/main/MarkovStateModel/MSM_validate.png?raw=true" alt="alt text" width="400">


Requirements: sda input file, trajectory file with related fcomplexes file.

### Selection of number of clusters and lagtime

In order to determine the number of cluster and lagtime to use for building the MSM, use the `Validate_MSM.py` python script which provides you a Within-Cluster Sum of Squares (WCSS) analysis together with a eigenvalues implied timescale analysis and a Chapman-Kolmogorov Test.

    python Validate_MSM.py --folder_prefix folder_trajectories --output_folder validate_folder --list_enc 3 6 16 46 48

### Extracting encountered trajectories from trajectory file

First you need to extract only the encountered trajectory from a trajectory file. The `Get_encounter_traj.py` python scrips extract the encounter trajectories by time reversing closest encounter complexes until com-com distance is over a user-defined cutoff. You can run it via:

    python  Get_encounter_traj.py sda.in ../data_grid/p1_noh.pdb ../data_grid/p2_noh.pdb 50

Output files are generated in the folder folder_ftrajectories

### Formatting encounered trajectories

The previous scripts extract the encountered trajectories in SDA format. First to make all the trajectories consistent between each other, you need to put all of them in the same reference system. To bring all the data in the reference system of the target protein use:

    python  Create_allign_enco_traj.py sda.in

Then it is also possible to convert the trajectories in xyz format to make them more accessible:

    python  Create_xyz_encoun_traj.py sda.in ../data_grid/p1_noh.pdb ../data_grid/p2_noh.pdb

### Build Markov State Model

Once you have the xyz trajectory files, it is possible to build the Markov State Models using `build_MSM_single_folder`. This python script makes use of `deeptime` python library

    python build_MSM_single_folder.py --folder_xyz folder_trajectories_1_xyz --folder_msm folder_MSM --num_clus 5 --max_iter_kmeans 50000

There are several other self-explaining flags that can be investigated via 

    python build_MSM_single_folder.py -h


## MSM from multiple trajectories

Most of the time SDA (especially in sdamm) doesn't have a single trajectory file where all the encountered trajectories are stored. If you have multiple encounter trajectory and complexes files which are not ordered, it is possible to use the bash script `encounter_traj_xyz_vmd.sh` which does the initial three steps for a list of encounter files. It is mandatory to specify in the bash script the list of encounter files that should be to analysed.

    ./encounter_traj_xyz_vmd.sh

the list of encounter trajectory files can be obtained via the `find_encounter_files.py` python script which can iterate over all complexes files and detect the ones which have a single or multiple encountered complexes. Run it providing the total number of soluter (include protein, ligand and crowders) and the prefix of the complexes file:

    python find_encounter_files.py 1107 sdamm_complexes_

The output is a list of integer which should be copied into `encounter_traj_xyz_vmd.sh`

Now, it is possible to build a MSM using all those input trajectories together. To do that run:

    python build_MSM.py --folder_xyz folder_trajectories --folder_msm folder_MSM --num_clus 5 --max_iter_kmeans 50000 --list_enc 84 159 190

where `list_enc` is again the list of integers provided by the `find_encounter_files.py` python script which specify the index of the encounter trajectories.


#### Other

The repository contains also a bash script to copy the encountered files (trajectory files and respective sda input and complexes files) into a new folder. This might be usefule if the simulations have run on a cluster, and the analysis should run on workstations: instead of remotly copying all the files, only the encounter files should be transfered. To do that, firts run `find_encounter_files.py` to get the list of indexes of the encountered files, copy the list into `copy_to_transfer.sh` file and run

    mkdir to_transfer
    ./copy_to_transfer sdamm_complexes_ trajectories_ sdamm_crowd_

This will copy the trajectory, complexes and input files with the provided prefixes into the folder `to_transfer`. Then it is possible to copy locally only that specific folder.
