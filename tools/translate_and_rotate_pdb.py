import numpy as np
from Bio import PDB
from scipy.spatial.transform import Rotation as R

def calculate_center_of_mass(coords, weights):
    weighted_coords = coords * weights[:, np.newaxis]
    com = weighted_coords.sum(axis=0) / weights.sum()
    return com

def translate_structure(structure, translation_vector):
    for atom in structure.get_atoms():
        atom.coord += translation_vector

def rotate_structure(structure, rotation_matrix):
    for atom in structure.get_atoms():
        atom.coord = np.dot(rotation_matrix, atom.coord)

def main(input_pdb, output_pdb, target_com, rotation_matrix=None, euler_angles=None):
    # Load the PDB structure
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    # Get all atom coordinates and weights (assumed to be atomic masses for simplicity)
    atom_coords = []
    atom_weights = []
    for atom in structure.get_atoms():
        atom_coords.append(atom.coord)
        atom_weights.append(atom.mass if atom.mass is not None else 1.0)  # Use atomic mass or default to 1.0

    atom_coords = np.array(atom_coords)
    atom_weights = np.array(atom_weights)

    # Calculate the current center of mass
    current_com = calculate_center_of_mass(atom_coords, atom_weights)

    # Compute the translation vector to move current COM to origin (for proper rotation)
    to_origin_vector = -current_com

    # Translate structure to move center of mass to the origin
    translate_structure(structure, to_origin_vector)

    # Rotate structure if a rotation matrix or Euler angles are provided
    if rotation_matrix is not None:
        rotate_structure(structure, rotation_matrix)
    elif euler_angles is not None:
        # Generate rotation matrix from Euler angles
        rot = R.from_euler('xyz', euler_angles, degrees=True)  # Assuming xyz Euler angles
        rotate_structure(structure, rot.as_matrix())

    # Compute the translation vector to move the center of mass to the target
    translation_vector = np.array(target_com)

    # Translate structure to the target center of mass
    translate_structure(structure, translation_vector)

    # Save the new structure to a new PDB file
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Translate and rotate a PDB structure to a specified center of mass and orientation.")
    parser.add_argument('input_pdb', help="Input PDB file")
    parser.add_argument('output_pdb', help="Output PDB file")
    parser.add_argument('target_com', type=float, nargs=3, help="Target center of mass (x, y, z)")
    
    # Option for Euler angles
    parser.add_argument('--euler_angles', type=float, nargs=3, help="Euler angles for rotation (roll, pitch, yaw) in degrees", default=None)
    
    # Option for direct rotation matrix
    parser.add_argument('--rotation_matrix', type=float, nargs=9, help="3x3 Rotation matrix (9 values in row-major order)", default=None)

    args = parser.parse_args()

    # Process rotation matrix if provided as argument
    if args.rotation_matrix:
        rotation_matrix = np.array(args.rotation_matrix).reshape(3, 3)
    else:
        rotation_matrix = None

    # Call the main function
    main(args.input_pdb, args.output_pdb, args.target_com, rotation_matrix=rotation_matrix, euler_angles=args.euler_angles)
