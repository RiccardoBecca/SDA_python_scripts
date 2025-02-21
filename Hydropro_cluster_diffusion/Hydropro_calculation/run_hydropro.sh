#!/bin/bash


#echo "This script will run a python script Hydropro.py to generate input file (hydropro.dat) for HYDROpro software for diffusion coefficient calculations."

#HYDROPRO=/hits/basement/mcm/munizcam/Hydropro
HYDROPRO=/hits/basement/mcm/beccarro/Hydropro
#PYTHON_SCRIPT=/hits/fast/mcm/munizcam/HIV-protease/BD/$2/hydropro

input=$1

#===========================================

export HYDROPRO
#export PYTHON_SCRIPT

#===========================================

python $HYDROPRO/Hydropro.py $1 > ${input%.*}-prepare_input.out

if [ $1 == "p2_noh.pdb" ]; then
       sed -i 's/2.9,            !AER, radius of primary elements/1.2,            !AER, radius of primary elements/' hydropro.dat
fi

${HYDROPRO}/hydropro10-lnx.exe hydropro.dat $1 > ${input%.*}-hydropro.out

python $HYDROPRO/Hydropro.py ${input%.*}-res.txt > ${input%.*}-sda_values.dat
