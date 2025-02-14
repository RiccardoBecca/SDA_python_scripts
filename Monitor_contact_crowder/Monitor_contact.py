import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import argparse
from scipy.spatial.distance import cdist
from tqdm import tqdm


def compute_center_of_geometry(X):
    """
    Computes the center of geometry for a given set of atomic positions.
    
    Parameters:
        X (numpy.ndarray): Atomic positions of shape (1, N, 3), where N is the number of atoms.
    
    Returns:
        numpy.ndarray: The center of geometry as a (3,) array.
    """
    return X.mean(axis=1).squeeze()  # Compute mean over atoms (axis=1) and remove extra dimensions



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


def get_position_sda(line_splitted, molecule, com_molecule):
    """
    Computes atomic positions in sda reference system
    
    Parameters:
        line_splitted : line in a trajectory file splitted
        molecule : md trajectory file object of the molecule
        com_molecule : molecule center of mass
    
    Returns:
        x : position of the molecule in sda ref system
    """
    x=molecule.xyz*10 -com_molecule[0]
    #Save traslations coordinates
    tx, ty, tz = float(line_splitted[2]), float(line_splitted[3]), float(line_splitted[4])
    
    #Save rotation Matrix elements
    r1x, r1y, r1z = float(line_splitted[5]), float(line_splitted[6]), float(line_splitted[7]) 
    r2x, r2y, r2z = float(line_splitted[8]), float(line_splitted[8]), float(line_splitted[10]) 
    r3x, r3y, r3z = vectorial_product([r1x, r1y, r1z],[r2x, r2y, r2z])
    
    #Save the rotation Matrix elements in the rotation Matrix
    RM=np.zeros((3,3))
    RM[0,:]=r1x,r1y,r1z
    RM[1,:]=r2x,r2y,r2z
    RM[2,:]=r3x,r3y,r3z
    
    #Transpose cause Fortran saves by column,
    RM=RM.transpose()

    x=x@RM.T+np.array([tx, ty, tz ])

    return x


def main(args):
    """
    Main function
    """
    sda_input_file=args.sda_input_file
    path_molecule=args.pdb_molecule
    path_crowd=args.pdb_crowders
    ftrajectories=args.trajectory
    dist_coms=args.dist_coms
    group_type=args.mol_number
    minimum_dist=args.contact_dist

    #read sda input
    L, total_solutes = read_sda_input(sda_input_file)

    #load pdb of the molecule
    molecule=md.load(path_molecule)
    molecule_top=molecule.topology
    com_molecule=md.compute_center_of_mass(molecule)*10

    #load pdb of the crowder
    crowd=md.load(path_crowd)
    crowd_top=crowd.topology
    com_crowd=md.compute_center_of_mass(crowd)*10


    #get array to monitor the contact and time of the contacts
    crowd_residue_contacts=np.zeros(crowd_top.atom(-1).residue.index+1) #+1 because pdb starts enumerating at 1 and python at 0
    crowd_number=np.zeros(total_solutes)
    time_bounded=np.zeros(total_solutes)
    previous_bounded=np.zeros(total_solutes)
    times=[]

    #Read trajectory file
    with open(ftrajectories,"r") as fin:
        traj_lines=fin.readlines()[2:]
    
    for i,line in tqdm(enumerate(traj_lines), total=len(traj_lines)):

        line_splitted=line.split()
    
        if int(line_splitted[1])==group_type:

            molecule_pos = get_position_sda(line_splitted, molecule, com_molecule)
            molecule_com=compute_center_of_geometry(molecule_pos)

            #translate the molecule in the center of the box (easy way to avoid counting PBC)
            vt=L/2-molecule_com
            molecule_pos=molecule_pos+vt
            

        if int(line_splitted[1])>1:   
            crowder_pos = get_position_sda(line_splitted, crowd, com_crowd)

            #translate atoms position by vt because we centered the molecule in the center of the box
            crowder_pos = crowder_pos +vt

            #print(np.sqrt(((compute_center_of_geometry(molecule_pos)-compute_center_of_geometry(crowder_pos))**2).sum()))

            #check if the two coms are withing a bigger distance (to save computational time)
            if np.sqrt(((compute_center_of_geometry(molecule_pos)-compute_center_of_geometry(crowder_pos))**2).sum()) < dist_coms:
                dist_matrix=cdist(molecule_pos[0],crowder_pos[0])
                min_dist=np.min(dist_matrix)
                _, atom_crowd=np.unravel_index(np.argmin(cdist(molecule_pos[0],crowder_pos[0])), cdist(molecule_pos[0],crowder_pos[0]).shape)
                if min_dist<=minimum_dist:
                    atom_selected=crowd_top.atom(atom_crowd)
                    residue_selected=atom_selected.residue.index
                    crowd_residue_contacts[residue_selected]+=1 #remember this is translated by 1:python starts enumerating at 1, pdb at 0
                    crowd_number[int(line_splitted[1])-1]+=1

                    #check if previous time was bounded
                    #if previous_bounded[int(line.split()[2])-1]==1:
                    time_bounded[int(line_splitted[1])-1]+=1
                    previous_bounded[int(line_splitted[1])-1]=1
            
                if min_dist>minimum_dist and previous_bounded[int(line_splitted[1])-1]==1:
                    times.append(time_bounded[int(line_splitted[1])-1])
                    previous_bounded[int(line_splitted[1])-1]=0
                    time_bounded[int(line_splitted[1])-1]=0

            if np.sqrt(((compute_center_of_geometry(molecule_pos)-compute_center_of_geometry(crowder_pos))**2).sum()) >= dist_coms and previous_bounded[int(line_splitted[1])-1]==1:
                times.append(time_bounded[int(line_splitted[1])-1])
                previous_bounded[int(line_splitted[1])-1]=0
                time_bounded[int(line_splitted[1])-1]=0

    np.save(args.output_folder+"/"+"crowd_residue_contacts_"+args.trajectory, crowd_residue_contacts)
    np.save(args.output_folder+"/"+"bounded_times_"+args.trajectory, np.array(times))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Python script to monitor contacts between molecule and crowders from sdamm trajectories",
                                     epilog="Example usage:\n"
                                    "python Monitor_contact.py --sda_input_file sdamm_crowd_1.in --trajectory trajectories_1"
                                    " --pdb_molecule p2_noh.pdb --pdb_crowders crowder_noh.pdb --mol_number 2 "
                                    "--dist_coms 30.0 --contact_dist 4.5 --output_folder crowder_contacts")
    parser.add_argument("--sda_input_file", help = "sda input file name", type=str, required = True)
    parser.add_argument("--trajectory", help = "filename for the trajecotry", type=str, required = True)
    parser.add_argument("--pdb_molecule", help = "pdb for the first molecule", type=str, required = True)
    parser.add_argument("--pdb_crowders", help = "pdb for the second molecule", type=str, required = True)
    parser.add_argument("--mol_number", help = "number of molecule in trajectory file", type=int, required = True)
    parser.add_argument("--dist_coms", help = "necessary distance to check if there is a contact (to speed up)", type=float, required = True)
    parser.add_argument("--contact_dist", help = "distance to define a contact", type=float, required = True)
    parser.add_argument("--output_folder", help = "name of the output folder where to save data", type=str, required = True)


    args = parser.parse_args()
    

    main(args)