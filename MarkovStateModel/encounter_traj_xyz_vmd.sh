#!/bin/sh

for i in 84 159 190 195 242 252 288 316 339 371 450 477 518 537 544 632 638 655 879 888 893 905 908 930 953
do
	echo Analyze file sdamm_crowd_$i.in
	python Get_encounter_traj.py sdamm_crowd_$i.in ../data_grid/p1_noh.pdb ../data_grid/p2_noh.pdb 30
	python Create_allign_enco_traj.py sdamm_crowd_$i.in
	python Create_xyz_encoun_traj.py sdamm_crowd_$i.in ../data_grid/p1_noh.pdb ../data_grid/p2_noh.pdb
done
