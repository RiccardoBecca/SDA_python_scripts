#!/bin/sh

subname_enc_files="$1"
subname_traj_files="$2"
subname_input_files="$3"


mkdir to_transfer

for i in 84 159 190 195 242 252 288 316 339 371 450 477 518 537 544 632 638 655 879 888 893 905 908 930 953
do
        cp $subname_enc_files$i to_transfer
        cp $subname_traj_files$i to_transfer
        cp $subname_input_files$i.in to_transfer
done
