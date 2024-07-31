### Tools for creating Markov State model from SDA.
Requirements: sda input file, trajectory file with related fcomplexes file.

#### Example usage

1. From trajectory extract encounter trajs by time reversing until com-com distance is over a user-defined cutoff.
    Input: sda input file (trajectory file and encounter complexes files should also be there), p1.pdb, p2.pdb and cutoff
   
    $ python  Get_encounter_traj.py sda.in ../data_grid/p1_noh.pdb ../data_grid/p2_noh.pdb 50
