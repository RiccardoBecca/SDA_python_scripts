import mdtraj as md
from tqdm import tqdm
import numpy as np
import sys
import os
import shutil
import warnings
warnings.filterwarnings("ignore")

def create_pdb_from_coordinates(coordinates, output_file):
    """
    Create a PDB file from a list of coordinates.
    
    Parameters:
    - coordinates: A list of tuples, where each tuple contains the (x, y, z) positions of an atom.
    - output_file: The name of the output PDB file.
    """
    
    with open(output_file, 'w') as pdb_file:
        for i, (x, y, z) in enumerate(coordinates, start=1):
            # Construct each PDB line according to the format
            pdb_line = (
                f"ATOM  {i:5d} CL   BEN A{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          CL  \n"
            )
            pdb_file.write(pdb_line)
            

def printUsage():

    print ("""
    NAME
        Create_xyz_encoun_traj.py

    DESCRIPTION

        Python script to write xyz file from extacted trajectory from Get_encounter_traj.py script

    ARGUMENTS

        It takes 3 parameters:
        1.)     Name of sda input file
        2.)     File/path of the p1_noh.pdb
        3.)     File/path of the p1_noh.pdb
        

    EXAMPLE

        python  Create_xyz_encoun_traj.py sda.in ../data_grid/p1_noh.pdb ../data_grid/p2_noh.pdb

    OUTPUT

        Output files with xyz coordinates of the solute 2 in the reference system of the p1 com are generated in folder_ftrajectories_xyz
        Output pdb files with xyz coordinates of the solute 2 in the reference system of the p1 com are generated in folder_ftrajectories_vmd: fast way to visualize trajs in vmd

    """)



if "-h" in sys.argv:
    print(printUsage())
    sys.exit()

with open(sys.argv[1], "r") as fin:
    lines_sda_input=fin.readlines()
for line in lines_sda_input:
    if "xmin" in line:
        xmin=float(line.split()[-1])
    if "xmax" in line:
        xmax=float(line.split()[-1])
    if "ftrajectories" in line:
        ftrajectories=line.split()[-1]

file_num=len(os.listdir(f"./folder_{ftrajectories}"))

L=xmax-xmin
D=L/2

#p2=md.load("p2_noh.pdb")
#p1=md.load("p1_noh.pdb")
p1_noh_file=sys.argv[2]
p2_noh_file=sys.argv[3]
p2=md.load(p2_noh_file)
p1=md.load(p1_noh_file)
p2_com=md.compute_center_of_geometry(p2)*10
p1_com=md.compute_center_of_geometry(p1)*10

if os.path.exists(f"folder_{ftrajectories}_xyz"):
    shutil.rmtree(f"folder_{ftrajectories}_xyz")
os.mkdir(f"folder_{ftrajectories}_xyz")


for t in tqdm(range(1,file_num)):
    subfile_traj=f"./folder_{ftrajectories}/{ftrajectories}_{t}"
    #subfile_traj=f"trajectories_1_{t}"
    with open(subfile_traj,"r") as fin:
        lines_fin=fin.readlines()


    with open(f"folder_{ftrajectories}_xyz/{ftrajectories}_{t}_xyz","w+") as fout:
    #with open(subfile_traj+"_xyz","w+") as fout:
        for line_fin in lines_fin[2:]:
            if int(line_fin.split()[2])==1:
                tx=float(line_fin.split()[3])
                ty=float(line_fin.split()[4])
                tz=float(line_fin.split()[5])
                r1x=float(line_fin.split()[6])
                r1y=float(line_fin.split()[7])
                r1z=float(line_fin.split()[8])
                r2x=float(line_fin.split()[9])
                r2y=float(line_fin.split()[10])
                r2z=float(line_fin.split()[11])
                r3x=r1y*r2z-r1z*r2y
                r3y=r1z*r2x-r1x*r2z
                r3z=r1x*r2y-r1y*r2x

                RP=np.zeros((3,3))
                RP[0,:]=r1x,r1y,r1z
                RP[1,:]=r2x,r2y,r2z
                RP[2,:]=r3x,r3y,r3z
                RP=RP.transpose()#magari da cancellare
                

            if int(line_fin.split()[2])==2:
                delta_tx=float(line_fin.split()[3])-tx
                delta_ty=float(line_fin.split()[4])-ty
                delta_tz=float(line_fin.split()[5])-tz


                delta_vec=np.array([delta_tx, delta_ty, delta_tz])
                delta_vec = delta_vec-L*np.round(delta_vec/L)

                r1xl=float(line_fin.split()[6])
                r1yl=float(line_fin.split()[7])
                r1zl=float(line_fin.split()[8])
                r2xl=float(line_fin.split()[9])
                r2yl=float(line_fin.split()[10])
                r2zl=float(line_fin.split()[11])
                r3xl=r1yl*r2zl-r1zl*r2yl
                r3yl=r1zl*r2xl-r1xl*r2zl
                r3zl=r1xl*r2yl-r1yl*r2xl
                Rl=np.zeros((3,3))
                Rl[0,:]=r1xl,r1yl,r1zl
                Rl[1,:]=r2xl,r2yl,r2zl
                Rl[2,:]=r3xl,r3yl,r3zl
                Rl=Rl.transpose()#magari da cancellare
                
                new_lig_pos=np.zeros((p2.xyz[0]).shape)
                
                for i,atom in enumerate(p2.xyz[0]):
                    new_lig_pos[i]=RP.transpose()@(Rl@(atom*10-p2_com[0])+delta_vec)
                    
                
                new_com_lig=np.zeros((3))
                for atom in new_lig_pos:
                    new_com_lig+=atom
                new_com_lig=new_com_lig/len(new_lig_pos)
                
                
                xsim='{0:.3f}'.format(new_com_lig[0]).rjust(7)
                ysim='{0:.3f}'.format(new_com_lig[1]).rjust(7)
                zsim='{0:.3f}'.format(new_com_lig[2]).rjust(7)
                
                newline=xsim+" "+ysim+" "+zsim+"\n"

                fout.write(newline)


if os.path.exists(f"folder_{ftrajectories}_vmd"):
    shutil.rmtree(f"folder_{ftrajectories}_vmd")
os.mkdir(f"folder_{ftrajectories}_vmd")

for t in tqdm(range(1,file_num)):
    coordinates = []
    with open(f"folder_{ftrajectories}_xyz/{ftrajectories}_{t}_xyz","r") as fin:
        lines=fin.readlines()
        for line in lines:
            coordinates.append((float(line.split()[0]), float(line.split()[1]), float(line.split()[2])))
    create_pdb_from_coordinates(coordinates, f"./folder_{ftrajectories}_vmd/output_1_{t}.pdb")

