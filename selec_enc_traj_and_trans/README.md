## Selection of encounter trajectories and copy them into a folder

Usually simulations are on the cluster. Then you want to download them to analyze or use MSM or whatever.
In this folder there are scripts to select the encounter trajectories from SDA and copy them into a folder.
Otherwise you will copy all the tarjectories but most of them will not be encountered trajs

### Example usage

  1. First step is to identify the encounter trajectories. You go through all the sdamm_complexes and find which are the files with some encountered complexes i.e. which are the files with not onl the starting position. You can use **find_encounter_files.py** fo tha. You should provide the total numer of solutes, the initial part of the encounter complexes file names, all te encounter complexes files.

      $ python find_encounter_files.py 2345 sdamm_complexes_ sdamm_complexes_{1..1000}

      This programm will give as output a string with all the files "number files" whichh contain a/some encounter complexes inside. You need to copy that string into the bash script (follow Step2).

  2. Once you have the string which tells which are the encounter files with a/some encounter complexes inside, you copy that string into the bash script **copy_to_transfer.sh** file. This bash script will create a new directory named **to_transfer** and will copy all the "encounterd" encounter files, trarjectories and input files. Be carefule, you need that the sda input files should end with '.in', otherwise it does not work. You should provide the bash script with the initial part of the names of the encounter_complexes. trajectories and input files.

       $ ./copy_to_transfer.sh sdamm_complexes_ trajectories_ sdamm_crowd_

     Now you have all the files useful for the MSM in the to_transfer file and you can scp only that folder.
