### Tools for creating Markov State model from SDA.
Requirements: sda input file, trajectory file with related fcomplexes file.

#### Example usage

1. From trajectory extract encounter trajs by time reversing until com-com distance is over a user-defined cutoff.
    Input: sda input file (trajectory file and encounter complexes files should also be there), p1.pdb, p2.pdb and cutoff
   
    $ python  Get_encounter_traj.py sda.in ../data_grid/p1_noh.pdb ../data_grid/p2_noh.pdb 50

    Output files are generated in the folder folder_ftrajectories

2. Python script to write rewrite encounter traj file from sda into new trajs files in sda format in the system of reference of p1_com
    Input: sda input file, also the extract trajectory from step 1 should be there.

   $ python  Create_allign_enco_traj.py sda.in

3. Python script to write xyz file from extracted trajectory from Get_encounter_traj.py script

   $ python  Create_xyz_encoun_traj.py sda.in ../data_grid/p1_noh.pdb ../data_grid/p2_noh.pdb

4. Build Markov State Model

    $ python build_MSM.py [-h] --folder_xyz FOLDER_XYZ --folder_msm FOLDER_MSM --num_clus NUM_CLUS [--lagtime LAGTIME] [--seed_kmeans SEED_KMEANS] [--max_iter_kmeans MAX_ITER_KMEANS] [--allow_disconnected ALLOW_DISCONNECTED] [--reversible REVERSIBLE]    
