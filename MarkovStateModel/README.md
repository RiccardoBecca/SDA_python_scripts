# Markov State Model for SDA.



<img src="https://github.com/RiccardoBecca/SDA_python_scripts/blob/main/MarkovStateModel/MSM_trypsin-ben.png?raw=true" alt="alt text" width="400">


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


#### Other

Most of the time SDA (especially in sdamm) doesn't have a single trajectory file where all the encountered trajectories are stored. If you have multiple encounter trajectory and complexes files which are not ordered, it is possible to use the bash script `encounter_traj_xyz_vmd.sh` which does the initial three steps for a list of encounter files. It is mandatory to specify in the bash script the list of encounter files that should be to analysed.

    ./encounter_traj_xyz_vmd.sh

You just need to change the list of the foor loop providing the correct list. The list can be obtained using the **find_encounter_files.py** in the other folder.

Of course now you should provide all the folders to build the MSM. The program **build_MSM.py** builds it and it needs the prefix name of the folder (which should end with "_xyz"), the name for the msm folder where outputs will be saved, the number of clusters to use and a list of integers wich indicates the name of the folders where to take the trajectories. For example

     python build_MSM.py --folder_xyz folder_trajectories --folder_msm folder_MSM --num_clus 5 --max_iter_kmeans 50000 --list_enc 84 159 190

will look into the folders folder_trajectories_84_xyz, folder_trajectories_159_xyz and folder_trajectories_190_xyz and load the data trajs there to build the MSM. This is intended to be the script where you have multiple folders with data. If you have a single folder with the data you can also use build_MSM_single_folder.py
