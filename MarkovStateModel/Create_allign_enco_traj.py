import shutil
import os
import sys
import numpy as np


def printUsage():

    print ("""
    NAME
        Create_allign_enco_traj.py

    DESCRIPTION

        Python script to write rewrite traj file from sda into new trajs files in sda format in the system of reference of p1_com

    ARGUMENTS

        It takes 1 parameters:
        1.)     Name of sda input file
        

    EXAMPLE

        python  Create_allign_enco_traj.py sda.in

    OUTPUT

        Output traj files in sda output format translated and rotated to have solutes in p1_com senter of geometry are written in folder_ftrajectories_center_sda format

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

L=xmax-xmin
D=L/2

file_num=len(os.listdir(f"./folder_{ftrajectories}"))

if os.path.exists(f"folder_{ftrajectories}_center_sda"):
    shutil.rmtree(f"folder_{ftrajectories}_center_sda")
os.mkdir(f"folder_{ftrajectories}_center_sda")

for t in range(1,file_num):
    subfile_traj=f"./folder_{ftrajectories}/{ftrajectories}_{t}"
    with open(subfile_traj,"r") as fin:
        lines_fin=fin.readlines()


    with open(f"./folder_{ftrajectories}_center_sda/{ftrajectories}_{t}_center", "w+") as fout:
    #with open(subfile_traj+"_center","w+") as fout:
        fout.write(lines_fin[0])
        fout.write(lines_fin[1])
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
                #mult_R=R.transpose()@R
                alligned_R=RP.transpose()@RP

                r1x_base=(alligned_R)[0,0]
                r1y_base=(alligned_R)[0,1]
                r1z_base=(alligned_R)[0,2]
                r2x_base=(alligned_R)[1,0]
                r2y_base=(alligned_R)[1,1]
                r2z_base=(alligned_R)[1,2]


                newtx='{0:.3f}'.format(0).rjust(7)
                newty='{0:.3f}'.format(0).rjust(7)
                newtz='{0:.3f}'.format(0).rjust(7)
                newr1x='{0:.3f}'.format(r1x_base).rjust(6)
                newr1y='{0:.3f}'.format(r1y_base).rjust(6)
                newr1z='{0:.3f}'.format(r1z_base).rjust(6)
                newr2x='{0:.3f}'.format(r2x_base).rjust(6)
                newr2y='{0:.3f}'.format(r2y_base).rjust(6)
                newr2z='{0:.3f}'.format(r2z_base).rjust(6)


                new_line_prot1=line_fin[:26]+newtx+line_fin[33:35]+newty+line_fin[42:44]+newtz+line_fin[51:54]+newr1x+line_fin[60:63]+newr1y+line_fin[69:72]+newr1z+line_fin[78:81]+newr2x+line_fin[87:90]+newr2y+line_fin[96:99]+newr2z+line_fin[105:]
                fout.write(new_line_prot1)
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


                #R_inv=R.transpose()

                delta_vec_rotated=RP.transpose()@delta_vec

                delta_tx_rotated=delta_vec_rotated[0]
                delta_ty_rotated=delta_vec_rotated[1]
                delta_tz_rotated=delta_vec_rotated[2]




                newtxl='{0:.3f}'.format(delta_tx_rotated).rjust(7)
                newtyl='{0:.3f}'.format(delta_ty_rotated).rjust(7)
                newtzl='{0:.3f}'.format(delta_tz_rotated).rjust(7)


                newR=RP.transpose()@Rl
                newR=newR.transpose()#magari da cancellare
                newr1xl=newR[0,0]
                newr1yl=newR[0,1]
                newr1zl=newR[0,2]
                newr2xl=newR[1,0]
                newr2yl=newR[1,1]
                newr2zl=newR[1,2]

                Rl1x='{0:.3f}'.format(newr1xl).rjust(6)
                Rl1y='{0:.3f}'.format(newr1yl).rjust(6)
                Rl1z='{0:.3f}'.format(newr1zl).rjust(6)
                Rl2x='{0:.3f}'.format(newr2xl).rjust(6)
                Rl2y='{0:.3f}'.format(newr2yl).rjust(6)
                Rl2z='{0:.3f}'.format(newr2zl).rjust(6)


                new_line_prot2=line_fin[:26]+newtxl+line_fin[33:35]+newtyl+line_fin[42:44]+newtzl+line_fin[51:54]+Rl1x+line_fin[60:63]+Rl1y+line_fin[69:72]+Rl1z+line_fin[78:81]+Rl2x+line_fin[87:90]+Rl2y+line_fin[96:99]+Rl2z+line_fin[105:]
                fout.write(new_line_prot2)





