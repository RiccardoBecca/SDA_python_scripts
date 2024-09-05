import sys
import numpy as np
from tqdm import tqdm
import mdtraj as md
import os
import shutil
import warnings
warnings.filterwarnings("ignore")

def dist_com2(v1,v2, L):
    d = abs(np.array(v1)-np.array(v2))
    D=L/2
    d = np.where(d<D,d,L-d)
    d= np.sqrt(np.sum(d**2,axis=-1))
    return d
def dist_com(v1,v2, L):
    delta_vec=v1-v2
    delta_vec = delta_vec-L*np.round(delta_vec/L)
    d= np.sqrt(np.sum(delta_vec**2,axis=-1))
    return d
    

def printUsage():

    print ("""
    NAME
        Get_encounter_traj.py

    DESCRIPTION

        Python script to extract the encounter trajectories from the trajectory file. Save only frames where com-com distance within a certain cutoff

    ARGUMENTS

        It takes 3 parameters:
        1.)     Name of sda input file
        2.)     File/path of the p1_noh.pdb
        3.)     File/path of the p1_noh.pdb
        4.)     cutoff value
        

    EXAMPLE

        python  Get_encounter_traj.py sda.in ../data_grid/p1_noh.pdb ../data_grid/p2_noh.pdb 50

    OUTPUT

        Output files are generated in the folder folder_ftrajectories

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
    if "fcomplexes" in line:
        fcomplexes=line.split()[-1]
    if "ftrajectories" in line:
        ftrajectories=line.split()[-1]
        
# Pass by argument
p1_noh_file=sys.argv[2]
p2_noh_file=sys.argv[3]
cut=float(sys.argv[4])
    
    
if os.path.exists(f"folder_{ftrajectories}"):
    shutil.rmtree(f"folder_{ftrajectories}")
os.mkdir(f"folder_{ftrajectories}")

    
L=xmax-xmin
D=L/2


# In[1]:


with open(fcomplexes,"r") as fin:
    complexes_lines=fin.readlines()[2:]
    
#We do not want the first traj also if it is encountered: might be 
#biased because the strarting point for first traj is not b-surf
prev_sim=0 # BE CAREFUL

#This is a dictionary with the encount traj num as keys (as it is display in trajectory file)
#and the first frame which appears in the sdamm_complexes file as argument
trajs_frames={}
complexed_linesp={}
complexed_linesl={}

for l,line in enumerate(complexes_lines):
    #if the num of the sim in the line is equal to the previous one
    #it means that the current frame is a continuation of the previous one
    #but we have already found the bound state so we don't care
    if int(line.split()[0]) != prev_sim and int(line.split()[1])!=0: #second if to avoi initial config

        #New encounter traj. Save that traj and this frame which is the first frame
        #if the traj close to the active state
        prev_sim=int(line.split()[0])
        trajs_frames[line.split()[0]]=int(line.split()[1])
        complexed_linesp[line.split()[0]]=complexes_lines[l]
        complexed_linesl[line.split()[0]]=complexes_lines[l+1]


# In[2]:


### TO DO ###
# p1_com to subtract


p2=md.load(p2_noh_file)
p1=md.load(p1_noh_file)
p2_com=md.compute_center_of_geometry(p2)*10
p1_com=md.compute_center_of_geometry(p1)*10


# In[3]:



#tot_trajs_num=1000

### CORRECT

# Now we create len()

with open(ftrajectories,"r") as fin_t:
    traj_lines=fin_t.readlines()

prev_traj=0
file_num=1
encountered_prot=False
print_traj=False
first_line=True
tempo_lines=[]

x_p,y_p,z_p=0,0,0
x_l,y_l,z_l=0,0,0

for l,line_t in enumerate(tqdm(traj_lines[2:])):
    
    #Per ottimizzare: metti qui if condition che guarda se la current traj è o meno in complexes, se no non ha senso
    if line_t.split()[0] in trajs_frames:
        #check if it is not a crowder
        if int(line_t.split()[2])<=2:
        

            #check if it is solute 1 (protein), in that case, save coordinates as we need them 
            #to compute COM with ligand and decide whether to save or not the traj
            if int(line_t.split()[2])==1:
                x_p,y_p,z_p=float(line_t.split()[3]), float(line_t.split()[4]), float(line_t.split()[5])
                x_l,y_l,z_l=float(traj_lines[2:][l+1].split()[3]), float(traj_lines[2:][l+1].split()[4]), float(traj_lines[2:][l+1].split()[5])

            #check if the current raw traj is the same as before, otherwise re-initialize 
            #the list to save the raws and updtae the number of trajs
            if int(line_t.split()[0])>prev_traj:
                tempo_lines=[]
                print_traj=False 
                prev_traj=int(line_t.split()[0])


            # We want to save the raw if the trajectory is in the complex file (trajs_frames)
            # and if we are not yet in the first encountered frame of that traj (trajs_frames[line_t.split()[0]])
            if line_t.split()[0] in trajs_frames and int(line_t.split()[1]) <=trajs_frames[line_t.split()[0]]:
                encountered_prot=True
                tempo_lines.append(line_t)

            #If the distance between protein and ligand is over a certain cutoff, we initialize again
            #Cause we want to discard those frames where the ligand has bound the protein to the 
            # "wrong" place and after that is gone away before finding the correct pocket.
            #if dist_com([x_l,y_l,z_l],[x_p,y_p,z_p],L)>70 and encountered_prot==True and print_traj==False:
            if dist_com(np.array([x_l,y_l,z_l]),np.array([x_p,y_p,z_p]),L)>cut:
                tempo_lines=[]
                first_line=True
                encountered_prot=False


            # If the current traj is in the encounter complexes file (trajs_frames), and if the currenct frame
            # is higher then the first encounter frame in the complexes file (trajs_frames[line_t.split())
            # it means that we have collected in tempo_lines all the frames of the current traj where the two
            # coms are under the cutoff (re-initialize if discontinuity) untill the encounter complexes.
            # Therefre it is time so save all of them.
            if line_t.split()[0] in trajs_frames and int(line_t.split()[1]) >trajs_frames[line_t.split()[0]]:  
                if print_traj==False:
                    print_traj=True
                    with open(f"./folder_{ftrajectories}/{ftrajectories}_{file_num}","w+") as fout:
                        # Save the first two standard lines
                        fout.write(traj_lines[0])
                        fout.write(traj_lines[1])
                        for tempo_line in tempo_lines:
                            fout.write(tempo_line)
                        # Save the first encountered complex in the traj
                        fout.write(complexed_linesp[line_t.split()[0]])
                        fout.write(complexed_linesl[line_t.split()[0]])
                    file_num+=1
                    first_line=True
                    line_t=[]



