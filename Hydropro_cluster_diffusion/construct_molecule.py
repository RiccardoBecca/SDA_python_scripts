import numpy as np
from Bio import PDB
import random
import os
import shutil
import argparse
from Bio.PDB import Chain
import string

def get_available_chain_id(used_ids):
    """Return a new chain ID that is not in used_ids."""
    for letter in string.ascii_uppercase:
        if letter not in used_ids:
            return letter
    raise ValueError("Ran out of chain IDs!")  # Unlikely unless >26 chains

def rename_chains(structure, used_ids):
    """Rename chains in structure to avoid conflicts."""
    for chain in structure.get_chains():
        if chain.id in used_ids:  # If chain ID is already in use
            new_id = get_available_chain_id(used_ids)
            used_ids.add(new_id)  # Mark this ID as used
            chain.id = new_id
        else:
            used_ids.add(chain.id)  # Track existing chain ID

def compute_center_of_geometry(structure):
    coords = [atom.coord for atom in structure.get_atoms() if atom.element != 'H']
    return np.mean(coords, axis=0)

def random_unit_vector():
    phi = random.uniform(0, 2 * np.pi)
    costheta = random.uniform(-1, 1)
    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])

def has_clash(structure1, structure2, threshold=1.5):
    atoms1 = [atom for atom in structure1.get_atoms() if atom.element != 'H']
    atoms2 = [atom for atom in structure2.get_atoms() if atom.element != 'H']
    for atom1 in atoms1:
        for atom2 in atoms2:
            if np.linalg.norm(atom1.coord - atom2.coord) < threshold:
                return True
    return False

def move_molecule(structure, translation):
    for atom in structure.get_atoms():
        atom.coord += translation
        

def combine_pdbs(mol1_pdb, mol2_pdb, output_pdb, step=2.0):
    parser = PDB.PDBParser(QUIET=True)
    mol1 = parser.get_structure("mol1", mol1_pdb)
    mol2 = parser.get_structure("mol2", mol2_pdb)

    center1 = compute_center_of_geometry(mol1)
    center2 = compute_center_of_geometry(mol2)

    direction = random_unit_vector()
    translation = direction * step
    move_molecule(mol2, translation - center2 + center1)

    while has_clash(mol1, mol2):
        translation += direction * step
        move_molecule(mol2, direction * step)

    io = PDB.PDBIO()
    combined_structure = PDB.Structure.Structure("combined")
    model = PDB.Model.Model(0)
    chain = PDB.Chain.Chain("A")  # Use a single chain for simplicity

    # Track max residue ID to avoid conflicts
    max_res_id = max(res.id[1] for res in mol1.get_residues())

    # Add all residues from mol1
    for residue in mol1.get_residues():
        new_residue = PDB.Residue.Residue(residue.id, residue.resname, residue.segid)
        for atom in residue.get_atoms():
            new_residue.add(atom)
        chain.add(new_residue)

    # Add all residues from mol2, assigning new residue IDs
    for residue in mol2.get_residues():
        max_res_id += 1  # Increment ID to ensure uniqueness
        new_res_id = (' ', max_res_id, ' ')  # Assign a new unique residue ID
        new_residue = PDB.Residue.Residue(new_res_id, residue.resname, residue.segid)
        for atom in residue.get_atoms():
            new_residue.add(atom)
        chain.add(new_residue)

    model.add(chain)
    combined_structure.add(model)

    io.set_structure(combined_structure)
    io.save(output_pdb)




def combine_pdbs_old(mol1_pdb, mol2_pdb, output_pdb, step=2.0):
    parser = PDB.PDBParser(QUIET=True)
    mol1 = parser.get_structure("mol1", mol1_pdb)
    mol2 = parser.get_structure("mol2", mol2_pdb)
    
    center1 = compute_center_of_geometry(mol1)
    center2 = compute_center_of_geometry(mol2)
    
    direction = random_unit_vector()
    
    translation = direction * step
    move_molecule(mol2, translation - center2 + center1)
    
    while has_clash(mol1, mol2):
        translation += direction * step
        move_molecule(mol2, direction * step)
    
    io = PDB.PDBIO()
    combined_structure = PDB.Structure.Structure("combined")
    model = PDB.Model.Model(0)  # Create a model to hold both molecules
    for chain in mol1.get_chains():
        model.add(chain)
    for chain in mol2.get_chains():
        model.add(chain)
    combined_structure.add(model)

    
    io.set_structure(combined_structure)
    io.save(output_pdb)
    
# Example usage:

def main(args):

    for i in range(1,args.max_mol_bounded+1):
        path=f"./{i}_bounded"
        
        #create folder where to save bounded molecules
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        
        if i==1:
            for n in range(args.num_replica):
                combine_pdbs(args.first_molecule, args.second_molecule, path+f"/{n+1}_rep.pdb", step=1.0)
        if i>1:
            for n,prev_molecule in enumerate(os.listdir(f"./{i-1}_bounded")):
                combine_pdbs(args.second_molecule, f"./{i-1}_bounded/"+prev_molecule, path+f"/{n+1}_rep.pdb", step=1.0)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Python script to build (random) bounded complexes between two given molecules",
                                     epilog="Example usage:\n"
                                     "python construct_molecule.py --first_molecule p2_noh.pdb --second_molecule crowder_noh.pdb --max_mol_bounded 4 --num_replica 5")
    parser.add_argument("--first_molecule", help = "pdb for the first molecule", type=str, required = True)
    parser.add_argument("--second_molecule", help = "pdb for the second molecule", type=str, required = True)
    parser.add_argument("--max_mol_bounded", help = "max number of molecule 2 bounded to molecule 1", type=int, required = True)
    parser.add_argument("--num_replica", help = "number of replica for building bounded complex", type=int, required = True)
    
    args = parser.parse_args()
    
    main(args)



