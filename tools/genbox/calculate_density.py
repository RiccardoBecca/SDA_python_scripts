import argparse

def calculate_density(box_size, molecular_weight, num_molecules):
    """
    Calculate the density (g/L) in a cubic box given the box size, molecular weight, and number of molecules.

    Parameters:
    box_size (float): Box side length in Angstroms (Å).
    molecular_weight (float): Molecular weight in Daltons (Da).
    num_molecules (int): Number of molecules in the box.

    Returns:
    float: Density in g/L.
    """
    # Constants
    avogadro_number = 6.02214076e23  # Avogadro's number, molecules/mol
    dalton_to_grams = 1.66053906660e-24  # 1 Dalton = 1.66053906660e-24 grams
    liter_to_cubic_angstrom = 1e27       # 1 Liter = 1e27 cubic Angstroms

    # Calculate volume of the box in liters
    box_volume_liters = (box_size**3) / liter_to_cubic_angstrom

    # Convert molecular weight to grams
    molecular_weight_grams = molecular_weight * dalton_to_grams

    # Calculate the total mass of the molecules in the box
    total_mass_grams = num_molecules * molecular_weight_grams

    # Calculate the density
    density = total_mass_grams / box_volume_liters

    return density

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate the density in a cubic box given the box size, molecular weight, and number of molecules.")
    parser.add_argument("-b", "--box_size", type=float, required=True, help="Box side length in Angstroms (Å)")
    parser.add_argument("-m", "--molecular_weight", type=float, required=True, help="Molecular weight in Daltons (Da)")
    parser.add_argument("-n", "--num_molecules", type=int, required=True, help="Number of molecules in the box")

    # Parse arguments
    args = parser.parse_args()

    # Perform calculation
    density = calculate_density(args.box_size, args.molecular_weight, args.num_molecules)
    print(f"Density: {density:.6f} g/L")

