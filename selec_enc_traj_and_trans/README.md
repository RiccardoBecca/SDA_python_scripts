## Selection of encounter trajectories and copy them into a folder

Usually simulations are on the cluster. Then you want to download them to analyze or use MSM or whatever.
In this folder there are scripts to select the encounter trajectories from SDA and copy them into a folder.
Otherwise you will copy all the tarjectories but most of them will not be encountered trajs

### Example usage

  1. First step is to identify the encounter trajectories. You go through all the sdamm_complexes and find which are the files with some encountered complexes i.e. which are the files with not onl the starting position. You can use *find_encounter_files.py* fo tha. You should provide the total numer of solutes, the initial part of the encounter complexes file names, all te encounter complexes files.
