This is a script to monitor the contact in a sdamm trajectory. Specify number of solute

Be careful at line 146. Here it is >1 because we only have ligand and crowders (used on trajectories produced for computing effective diffusion coefficient). You want to convert into >2 in the case you have p1, p2 and crowder

Same for --mol_number. Here is 1 but for p1, p2 and crowder should be 2 if you want select the ligand
