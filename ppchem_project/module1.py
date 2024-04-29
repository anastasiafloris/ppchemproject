#calculating the age of a sample using the Carbon-14 isotope knowing that the half life of carbon-14 is 5700 years
    

 # Half-life of Carbon-14 in years
import math
half_life_c14 = 5700

    # Assuming remaining ratio of Carbon-14 in the sample
remaining_ratio = float(input("Enter the remaining ratio of Carbon-14 in the sample (between 0 and 1): "))

# Calculate the age using the formula: age = -(half_life) * ln(remaining_ratio) / ln(2)
age = -(half_life_c14) * math.log(remaining_ratio) / math.log(2)
print (age)

print(f"The age of the sample is approximately {age:.2f} years.")

# project idea: Create a Python package that allows users to analyze chemical compounds, perform molecular computations, visualize molecular structures, and generate various plots related to chemical properties.
# 1. Molecular Descriptors Calculation: RDKit to calculate molecular descriptors such as molecular weight, number of atoms, number of bonds, etc.
# 2. Chemical Similarity Search: use of functions to perform chemical similarity searches based on molecular fingerprints or other similarity metrics.
# 3. Molecule Drawing: Allows users to draw chemical structures using RDKit and visualize them using matplotlib.
# 4. Data Import and Export: Providing functions to import chemical data from various file formats (e.g., SDF, SMILES) and export results to CSV or other formats.
# 5. Visualization: Creation of functions to generate plots and visualizations of chemical properties, such as histograms of molecular weights, scatter plots of property distributions, or 3D molecular structures.

# use of RDKit to calculate molecular weight and molecular projection
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors

def main():
    # Prompt the user to enter a SMILES string
    smiles = input("Enter the SMILES string of the molecule: ")

    # Create a molecule from the user's SMILES input
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        # Calculate molecular weight
        mw = Descriptors.MolWt(mol)
        print("Molecular weight:", mw)

        # Draw the molecule
        img = Draw.MolToImage(mol)
        img.show()  # Show the image in an external viewer
    else:
        print("Invalid SMILES string")

if __name__ == "__main__":
    main()



    


