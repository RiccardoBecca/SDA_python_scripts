Python scripts for computing number of molecules required to target a density, and to compute number of molecules required to match specific density.

#### From box size (A side), molecular weight (Da), target density (g/L) to -> number of molecules

     python calculate_molecules_cli.py -d 75.0 -m 14000 -b 240

This tells you number of molecules required. Approximate the number to the next integer.


#### From box size (A side), molecular weight (Da), number of molecules to -> target density (g/L)

     python calculate_density.py -b 240 -m 720 -n 71

This is the other way around
