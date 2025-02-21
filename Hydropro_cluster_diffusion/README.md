## Translational diffusion constants

Thgis repository contains the code to perform back-of-the-envelope calculation (based on hydropro results) for how much the translational diffusion constant would change if the diluted molecule binds to a second molecule i.e. consider the bounded complex as a single bigger molecule

First, given the pdb files of the two molecules, the bounded complexes between the two should be created. This can be done using `construct_molecule.py` script. It takes randomly sample a direction and virtually moves from molecule 1 center of geometry in that direction with step of 1Ang. At each step, it tries to place the second molecule. It checks if there are clashes: if so it continuos moving along the direction, otherwise the bounded complex is retruned in a pdb file. Run:

     python construct_molecule.py --first_molecule mol1.pdb --second_molecule mol2.pdb --max_mol_bounded n_bounded --num_replica num_replica

here `n_bounded` is the maximum number of mol2 molecules bounded to mol1, and `num_replica` is the number of complexes generated for each number of bounded complexes.

The script has generated as many folder as the maximum molecule bounded requested, and for each number of molecule bounded, different replicas of the bounded complex with differently sample random directions. Now it is possible to compute translational diffusion constants for each of these complexes using Hydropro:

     ./Hydropro_calculation/run_all.sh

The script generates the result in the `diffusion_coefficients` folder. It is finally possible to plot the values using:

     python plot_diffusion_coefficients.py --max_mol_bounded 6 --diff_coeff_folder diffusion_coefficients/ --num_replica 10 --Dzero diff_zero

here diff_zero is the translational diffusion constant of mol1 on its own (i.e in water).
