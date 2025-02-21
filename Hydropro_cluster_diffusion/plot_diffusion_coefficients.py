import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

def main(args):

    diff_coeff=np.zeros((args.max_mol_bounded, args.num_replica))
    diff_coeff_means=np.zeros((args.max_mol_bounded+1))
    diff_coeff_std=np.zeros((args.max_mol_bounded+1))

    #Fill first value wit Dzero
    diff_coeff_means[0]=args.Dzero
    diff_coeff_std[0]=0.

    for n in range(args.max_mol_bounded):
        for f,file in enumerate(os.listdir(args.diff_coeff_folder+f"/{n+1}_bounded")):
            value=np.loadtxt(args.diff_coeff_folder+f"/{n+1}_bounded/"+file)
            diff_coeff[n,f]=value
        diff_coeff_means[n+1]=np.mean(diff_coeff[n,:])
        diff_coeff_std[n+1]=np.std(diff_coeff[n,:])

    plt.figure()
    plt.errorbar(x=np.arange((args.max_mol_bounded+1)), y=diff_coeff_means, yerr=diff_coeff_std)
    plt.scatter(np.arange((args.max_mol_bounded+1)), diff_coeff_means)
    plt.show()
        


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="plotting diffusion coefficients calculated with Hydropro from virtually clustered molecules")
    parser.add_argument('--max_mol_bounded', type=int, help="max number of molecule 2 bounded to molecule 1",required=True)
    parser.add_argument('--diff_coeff_folder', type=str, help="name of the folder where diffusion coefficients are saved",required=True)
    parser.add_argument('--num_replica', type=int, help="number of replica for building bounded complex",required=True)
    parser.add_argument('--Dzero', type=float, help="diffusion coefficient without any further bounded molecule",required=True)

    args=parser.parse_args()
    
    main(args)