import argparse

def calculate_molecules(density, molecular_weight, box_size):
    """
    Calculate the number of molecules required to achieve a specified density in a cubic box.

    Parameters:
    density (float): Density in g/L.
    molecular_weight (float): Molecular weight in Daltons (Da).
    box_size (float): Box side length in Angstroms (Å).

    Returns:
    int: Number of molecules needed.
    """
    # Constants
    avogadro_number = 6.02214076e23  # Avogadro's number, molecules/mol
    dalton_to_grams = 1.66053906660e-24  # 1 Dalton = 1.66053906660e-24 grams
    liter_to_cubic_angstrom = 1e27       # 1 Liter = 1e27 cubic Angstroms

    # Calculate volume of the box in liters
    box_volume_liters = (box_size**3) / liter_to_cubic_angstrom

    # Calculate the mass of the substance in the box (density × volume)
    total_mass_grams = density * box_volume_liters

    # Convert molecular weight to grams
    molecular_weight_grams = molecular_weight * dalton_to_grams

    # Calculate the number of molecules
    num_molecules = total_mass_grams / molecular_weight_grams

    return num_molecules

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate the number of molecules needed to match a given density in a cubic box.")
    parser.add_argument("-d", "--density", type=float, required=True, help="Density in g/L")
    parser.add_argument("-m", "--molecular_weight", type=float, required=True, help="Molecular weight in Daltons (Da)")
    parser.add_argument("-b", "--box_size", type=float, required=True, help="Box side length in Angstroms (Å)")

    # Parse arguments
    args = parser.parse_args()

    # Perform calculation
    num_molecules = calculate_molecules(args.density, args.molecular_weight, args.box_size)
    print(f"Number of molecules needed: {num_molecules}")

