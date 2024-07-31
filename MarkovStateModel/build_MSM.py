import sys
import numpy as np
import os
import shutil
from deeptime.markov.msm import MaximumLikelihoodMSM
from tqdm.notebook import tqdm  # progress bar (optional)
from deeptime.clustering import KMeans
import matplotlib.pyplot as plt
import argparse
from deeptime.plots import plot_markov_model
import matplotlib as mpl
import networkx as nx


def printUsage():

    print ("""
    NAME
        build_MSM.py

    DESCRIPTION

        Python script to build MSM from xyz files

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
    
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def create_pdb_from_coordinates(coordinates, output_file):
    """
    Create a PDB file from a list of coordinates.
    
    Parameters:
    - coordinates: A list of tuples, where each tuple contains the (x, y, z) positions of an atom.
    - output_file: The name of the output PDB file.
    """
    
    with open(output_file, 'w+') as pdb_file:
        for i, (x, y, z) in enumerate(coordinates, start=1):
            # Construct each PDB line according to the format
            pdb_line = (
                f"ATOM  {i:5d} CL   BEN A{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          CL  \n"
            )
            pdb_file.write(pdb_line)



parser = argparse.ArgumentParser(description = "Description for my parser")
parser.add_argument("--folder_xyz", help = "Folder which contains xyz files", type=str, required = True)
parser.add_argument("--folder_msm", help = "Folder where to save MarkovStateModel, it exists, overwrite", type=str,required = True)
parser.add_argument("--num_clus", help = "Number of clusters defined for KMeans evaluation", type=int,required = True)
parser.add_argument("--lagtime", help = "If provided, xyz trajs will be taken with the lagtime provided", type=int, required = False, default=1)
parser.add_argument("--seed_kmeans", help = "Set KMans seed", type=int, required = False, default=13)
parser.add_argument("--max_iter_kmeans", help = "Max iteration in KMeans cluster optimization.", type=int, required = False, default=500)
#parser.add_argument("--allow_disconnected", help = "If set, allow disconnected states in the MSM (not recommended)", action=argparse.BooleanOptionalAction)
parser.add_argument("--allow_disconnected", help = "If true, allow disconnected states in the MSM (not recommended)", type=str2bool, default=False)
parser.add_argument("--reversible", help = "If true, MSM constrained to satisfy detail balanced, default=False", type=str2bool, default=False)
                
argument = parser.parse_args()

folder_xyz=argument.folder_xyz
folder_msm=argument.folder_msm
num_clus=argument.num_clus
lagtime_msm=argument.lagtime
seed_kmeans=argument.seed_kmeans
max_iter_kmeans=argument.max_iter_kmeans
allow_disconnected_input=argument.allow_disconnected
reversible_input=argument.reversible

print()
print(f"Running KMeans with n_clusters = {num_clus}")
print(f"Running MSM with allow_disconnected = {allow_disconnected_input}")
print(f"Running MSM with reversibility = {reversible_input}")

data=[]
for file in os.listdir(folder_xyz):
    if "trajectories_1" in file:
        data.append(np.loadtxt(f"{folder_xyz}/{file}"))

estimator = KMeans(
    n_clusters=num_clus,  # place 100 cluster centers
    init_strategy='kmeans++',  # kmeans++ initialization strategy
    max_iter=0,  # don't actually perform the optimization, just place centers
    fixed_seed=seed_kmeans
)

clustering = estimator.fit(np.concatenate(data)).fetch_model()
assignments = clustering.transform(np.concatenate(data))
estimator.initial_centers = clustering.cluster_centers
estimator.max_iter = max_iter_kmeans
clustering_new = estimator.fit(np.concatenate(data)).fetch_model()

if os.path.exists(folder_msm):
    shutil.rmtree(folder_msm)
os.mkdir(folder_msm)


for i,cluster in enumerate(clustering_new.cluster_centers):
    coordinates=[(cluster[0], cluster[1], cluster[2])]
    create_pdb_from_coordinates(coordinates, f"./{folder_msm}/cluster_{i}.pdb")

    

plt.figure()
plt.loglog(clustering_new.inertias)
plt.xlabel("iteration")
plt.ylabel("inertia")
plt.title("KMeans++ inertia during training")
plt.savefig(f"{folder_msm}/inertia_kmean_clust")


dtrajs = [estimator.transform(t) for t in data]

msm = MaximumLikelihoodMSM(lagtime=lagtime_msm, reversible=reversible_input, allow_disconnected=allow_disconnected_input).fit_fetch(dtrajs)

if allow_disconnected_input==False and msm.n_states != num_clus:
    print("ERROR:")
    print("     Set allow_disconnected=False but found disconnected states")
    print("     Decrease number of cluster (Suggested)")
    print("     Set allow_disconnected=True (Not suggested)")
    sys.exit(2)
    
with open(f"./{folder_msm}/msm_rates", "w+") as fout:
    fout.write("Markov State Models transition matrix"+"\n")
    fout.write(f"Number of states = {msm.n_states}"+"\n")
    for i in range(msm.n_states):
        fout.write(f"state = {i}"+"\n")
        for j in range(msm.n_states):
            fout.write(f"p_{i}->{j} = {msm.transition_matrix[i,j]}"+"\n")
            
with open(f"./{folder_msm}/msm_mfpt", "w+") as fout:
    fout.write("Markov State Models mean first passage time matrix"+"\n")
    fout.write(f"Number of states = {msm.n_states}"+"\n")
    for i in range(msm.n_states):
        fout.write(f"state = {i}"+"\n")
        for j in range(msm.n_states):
            fout.write(f"mfpt_{i}->{j} = {msm.mfpt(i,j)}"+"\n")
            
plt.figure(figsize=(15,15))
plot_markov_model(msm)
plt.savefig(f"./{folder_msm}/msm_network")
