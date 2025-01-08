import os
import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse

import numpy as np

def unwrap_positions(positions, box_xmin=0, box_xmax=240, box_ymin=0, box_ymax=240, box_zmin=0, box_zmax=240):
    """
    Unwrap particle positions for a 3D box considering periodic boundary conditions,
    explicitly handling x, y, and z axes separately.
    
    Parameters:
        positions (numpy array): Array of shape (T, 3), where T is the number of timesteps.
        box_xmin, box_xmax: Bounds of the box in the x direction.
        box_ymin, box_ymax: Bounds of the box in the y direction.
        box_zmin, box_zmax: Bounds of the box in the z direction.
    
    Returns:
        numpy array: Unwrapped positions of shape (T, 3).
    """
    # Initialize unwrapped positions
    unwrapped_positions = np.copy(positions)
    T, d = positions.shape  # Number of timesteps (T) and dimensions (d=3 for x, y, z)
    
    # Box size in each direction
    box_size_x = box_xmax - box_xmin
    box_size_y = box_ymax - box_ymin
    box_size_z = box_zmax - box_zmin
    
    # Initialize crossing counters for x, y, z
    counter_x = 0
    counter_y = 0
    counter_z = 0
    
    for t in range(1, T):
        # Compute displacements for x, y, z
        dx = positions[t, 0] - positions[t - 1, 0]
        dy = positions[t, 1] - positions[t - 1, 1]
        dz = positions[t, 2] - positions[t - 1, 2]
        
        # Check for x boundary crossings
        if dx > box_size_x / 2:
            counter_x -= 1  # Crossed left to right
        elif dx < -box_size_x / 2:
            counter_x += 1  # Crossed right to left
        
        # Check for y boundary crossings
        if dy > box_size_y / 2:
            counter_y -= 1  # Crossed bottom to top
        elif dy < -box_size_y / 2:
            counter_y += 1  # Crossed top to bottom
        
        # Check for z boundary crossings
        if dz > box_size_z / 2:
            counter_z -= 1  # Crossed front to back
        elif dz < -box_size_z / 2:
            counter_z += 1  # Crossed back to front
        
        # Update unwrapped positions
        unwrapped_positions[t, 0] = positions[t, 0] + counter_x * box_size_x
        unwrapped_positions[t, 1] = positions[t, 1] + counter_y * box_size_y
        unwrapped_positions[t, 2] = positions[t, 2] + counter_z * box_size_z

    return unwrapped_positions


def read_sda_input_file(folder, file):
    """
    Read sda input file
    
    Parameters:
        - folder name selected
        - file name of a trajectory file
    
    Returns:
        timestep, frequency of printing and box size of the simulation
    """
    num_traj=file.replace('trajectories','')
    with open(f"./{folder}/sdamm_crowd"+num_traj+".in") as sda_input:
        sda_input_lines=sda_input.readlines()
    for sda_line in sda_input_lines:
        tempo_sda_line=sda_line.strip().split()
        if "dt1 = " in sda_line:
            dt=float(tempo_sda_line[2])
        if "freq_print = " in sda_line:
            frequency = int(tempo_sda_line[2])
        if "zmax = " in sda_line:
            box_size = float(tempo_sda_line[2])
            
    return dt, frequency, box_size
    
    
def load_positions(lines, solute_idx, start_time):
    """
    load position from trajectories. Only solute_idx positions are loaded
    
    Parameters:
        - lines of trjactory file
        - solute idx of interested solute
        - start time value for considering equilibration
    
    Returns:
        numpy array with positions
    """
    positions=[]
    
    # iterate over line of the trajectory file
    for line in lines[2:]:
        parts = line.split()  # Split the line only once
        current_solute_idx = int(parts[1])

        # if the solute index of that file is 1 i.e. the molecule of interest, save pos
        if current_solute_idx == solute_idx:
            pos_t = np.array([float(parts[2]), float(parts[3]), float(parts[4])])

            # save only positions above the start time
            if int(parts[0])>=start_time:
                positions.append(pos_t)
    positions=np.array(positions)
    
    return positions

def write_diffusion_coefficients(Ds, densities, name_output_file):
    """
    write diffusion coefficients into a file
    
    Parameters:
        - array with diffusion coefficients
        - array with density values
    
    """
    with open(name_output_file, "w+") as fout:
        for d,density in enumerate(densities):
            fout.write(str(density) + " " + str(Ds[d]) + "\n")

def main(rejected_frames, solute_idx, min_time_fit, max_time_fit, name_output_file, folder_figures):
    densities=[]
    Ds=[]

    # iterate over folders
    for folder in os.listdir("./"):
        #select only the folder with association simulations i.e. "assoc_" string in the name
        if "assoc_" in folder:
        #if "assoc_200gL" in folder or "assoc_100gL" in folder or "assoc_0gL" in folder:
            #count number of traj in that folder -> used as normalization constant in msd
            count_trajs=0
            #flag value to monitor if a first array should be defined
            first_traj=True
            
            # iterate over files in folder
            for file in tqdm(os.listdir(f"./{folder}")):
                #select only trajectory files
                if "trajectories" in file:
                    
                    # define position and unwrapped_position lists
                    positions = []
                    positions_unwrapped = []
                    with open(f"./{folder}/{file}", "r") as fin:
                        lines = fin.readlines()
                        
                    #Read parameters from sda_input file
                    dt, frequency, box_size = read_sda_input_file(folder,file)
                                       
                    # discard first start_time frame for equilibration
                    start_time=rejected_frames*frequency
                    time_step=dt*frequency
                    
                    #increment number of trajectories for that density
                    count_trajs+=1

                    # Load positions of requires solute
                    positions=load_positions(lines, solute_idx, start_time)

                    # Unwrap the trajectory
                    unwrapped_positions = unwrap_positions(positions, box_xmin=0, box_xmax=box_size, box_ymin=0, box_ymax=box_size, box_zmin=0, box_zmax=box_size)

                    # Number of timesteps and dimensions
                    T, d = unwrapped_positions.shape
                    
                    # if this is the first trajectory of that density, intialise the msd aray
                    if first_traj == True:
                        #remove 1 because frame 0 cannot count for msd
                        msd = np.zeros(T-1)
                        first_traj = False
                    
                    #compute msd of that traj
                    for t in range(1,T):
                        displacement=unwrapped_positions[t] - unwrapped_positions[0]
                        msd[t-1]+=np.sum(displacement**2)

            
            # average msd values over all trajectories
            times=[]
            for t in range(1,T):
                msd[t-1]/=count_trajs
                times.append(t*time_step)
            times=np.array(times)
            
            # set range of time (time unit) where to compute Diffusion coefficient
            long_time_regime = slice(min_time_fit, max_time_fit)
            
            # fit the msd vs time linear funciton
            slope, intercept = np.polyfit(times[long_time_regime], msd[long_time_regime], 1)

            #extract Diffusion coefficient value using Einstein relation
            D = slope / (2 * d)

            Ds.append(D)
            densities.append(int(folder.replace("assoc_","").replace("gL","")))
            print(f"Diffusion Coefficient {folder}: {D:.4f}")

    plt.figure()
    plt.scatter(np.array(densities),np.array(Ds), s=40)
    plt.plot(np.array(densities),np.array(Ds))
    plt.xlabel("Densities [g/L]")
    plt.ylabel("Diffusion coefficients [Ang^2/ps]")
    plt.savefig(folder_figures+"/diff_coeff_vs_dens")
    
    write_diffusion_coefficients(Ds, densities, name_output_file)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Comput diffusion coefficients for different densities")
    parser.add_argument("--rejected_frames", type=int, required=False, default=10, help="number of initial frames to reject to consider equilibration")
    parser.add_argument("--solute_idx", type=int, required=False, default=1, help="Solute index of which computing the Diffusion coefficient")
    parser.add_argument("--min_time_fit", type=int, required=False, default=0, help="minimum time from where to start msd vs t fit")
    parser.add_argument("--max_time_fit", type=int, required=False, default=-1, help="maximum time where to end msd vs t fit")
    parser.add_argument("--name_output_file", type=str, required=False, default="diff_coeff.txt", help="name of output file to save diffusion coefficients")
    parser.add_argument("--folder_figures", type=str, required=False, default="images_diff", help="folder where to save diffusion images")

    # Parse arguments
    args = parser.parse_args()

    # Perform calculation
    main(args.rejected_frames, args.solute_idx, args.min_time_fit, args.max_time_fit, args.name_output_file, args.folder_figures)
