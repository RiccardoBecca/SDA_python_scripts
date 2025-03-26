import argparse
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def read_sda_input(sda_input_file):
    """
    Reads sda input file
    
    Parameters:
        sda_input_file: name of sda input file
    
    Returns:
        L : box lenght
        total_solutes = number of solutes
    """
    with open(sda_input_file, "r") as fin:
        lines_sda_input=fin.readlines()
    for line in lines_sda_input:
        if "xmin" in line:
            xmin=float(line.split()[-1])
        if "xmax" in line:
            xmax=float(line.split()[-1])
        if "total_solutes" in line:
            total_solutes=int(line.split()[-1])

    return xmax-xmin, total_solutes

def distance_pbc(x0, x1, dimensions):
    """
    Computes PBC distances from relative distances
    
    Parameters:
        x0 and x1: two arrays. In this case x0.shape = (N_atoms, 3); x0.shape = (3)
        dimensions: box size
    
    Returns:
        numpy.ndarray: distances in PBC
    """
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    
    return np.sqrt((delta ** 2).sum(axis=-1))

def compute_center_of_geometry(X):
    """
    Computes the center of geometry for a given set of atomic positions.
    
    Parameters:
        X (numpy.ndarray): Atomic positions of shape (N, 3), where N is the number of atoms.
    
    Returns:
        numpy.ndarray: The center of geometry as a (3,) array.
    """

    return X.mean(axis=0).squeeze()  # Compute mean over atoms (axis=0) and remove extra dimensions

def read_frame(line):
    """
    Function to read SDA frame
    
    Parameters:
        line : SDA trajectory line
        
    Returns : frame, solute_idx, tx, ty, tz, r1x, r1y, r1z, r2x, r2y, r2z
        
    """
    parts=line.split()
    return int(parts[1]), int(parts[2]), float(parts[3]), float(parts[4]), float(parts[5]), float(parts[6]), float(parts[7]), float(parts[8]), float(parts[9]), float(parts[10]), float(parts[11])

def vectorial_product(vec1, vec2):
    """
    Computes vectorial product between two vectors
    
    Parameters:
        vec1 : list of coordinates of first vector
        vec2 : list of coordinates of second vector
    
    Returns:
        vec3 : vector defined by the vector product between vec1 and vec2
    """
    return vec1[1]*vec2[2]-vec1[2]*vec2[1],vec1[2]*vec2[0]-vec1[0]*vec2[2],vec1[0]*vec2[1]-vec1[1]*vec2[0]


def main(args):
    """
    Main function
    """

    #Read SDA input file
    L, total_solutes = read_sda_input(args.sda_input_file)

    #read reaction coordinate files
    with open(args.reaction_file, "r") as fin:
        lines= fin.readlines()
    
    #save reaction coordinate atoms in lists
    reaction_atoms_idx=[]
    indexes_saved=set()
    for line in lines:
        idx=int(line.split()[2])
        if idx not in indexes_saved:
            reaction_atoms_idx.append(idx)
            indexes_saved.add(idx)

    print(f"\nReaction atom indexes: {reaction_atoms_idx} \n")

    #load pqr file
    protein = mda.Universe(args.protein)
    protein_com=compute_center_of_geometry(protein.atoms.positions)

    # save crystallography positions of reaction atoms in list
    reaction_atoms_coordinates=[]
    for id in reaction_atoms_idx:
        print(protein.atoms[id -1]) #-1 is required because python is 0 indexed
        reaction_atoms_coordinates.append(protein.atoms[id -1].position)
    reaction_atoms_coordinates=np.array(reaction_atoms_coordinates)
        
    #Load trajectories
    with open(args.trajectory,"r") as fin:
        traj_lines=fin.readlines()[2:]
        
        
    min_distaces=[]
        
    #iterate over frames
    for line in tqdm(traj_lines):
        frame, solute_idx, tx, ty, tz, r1x, r1y, r1z, r2x, r2y, r2z = read_frame(line)
        r3x, r3y, r3z = vectorial_product([r1x, r1y, r1z],[r2x, r2y, r2z])
        
        #If the molecule is the protein
        if solute_idx == 1:
            #define rotation matrix
            RM = np.array([[r1x,r1y,r1z],[r2x, r2y, r2z], [r3x, r3y, r3z]])
            #Transpose cause Fortran saves by column,
            RM=RM.transpose()
            
            # x are coordinates of the entire protein, xrxna coordinates of reaction atoms
            x=np.dot(protein.atoms.positions-protein_com, RM.T)+np.array([tx, ty, tz])
            xrxna=np.dot(reaction_atoms_coordinates-protein_com, RM.T)+np.array([tx, ty, tz])
            
            distances_from_crowders=[]
            
        #else is a crowder
        else:
            #get crowder distances to reaction atoms coordinates
            distances=xrxna-np.array([tx, ty, tz])
            distances_pbc=distance_pbc(xrxna, np.array([tx, ty, tz]),L)
            distances_from_crowders.append(distances_pbc)
        
        # If we are at the last crowder, save the minimum distances between reactive atoms
        # and crowders
        if solute_idx == total_solutes:
            min_distaces.append(np.min(distances_from_crowders))

    plt.hist(min_distaces, bins=50)
    plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Python script to monitor contacts between C60 and residues involved in trypsin pocket",
                                     epilog="Example usage:\n"
                                    "python Monitor_distance.py --trajectory ./assoc_0.05/trajectories_1 --reaction_file ./assoc_0.05/p2_methyl.rxna --protein ./data/p1.pqr --sda_input_file ./assoc_0.05/sdamm_crowd_1.in"
                                    )
    parser.add_argument("--trajectory", help = "filename for the trajecotry", type=str, required = True)
    parser.add_argument("--reaction_file", help = "reaction coordinate file", type=str, required = True)
    parser.add_argument("--protein", help = "pqr file of protein", type=str, required = True)
    parser.add_argument("--sda_input_file", help = "sda input file name", type=str, required = True)

    args = parser.parse_args()
    

    main(args)
