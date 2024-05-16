# project idea: Create a Python package that allows users to analyze chemical compounds, perform molecular computations and visualize molecular structures.
# 1. Molecular Descriptors Calculation: RDKit to calculate molecular descriptors such as molecular weight, number of atoms and number of bonds.
# 2. Molecule Drawing: Allows users to draw chemical structures using RDKit.
# 3. Data Import and Export: Providing functions to import chemical data as CSV files. 
# 4. Visualization: 3D molecular interactive visulaization.
# 5. Presentation: displaying the molecular characteristics using pandas.

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

#the following code takes the 3D molecular structure of the molecule into account:
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt

def validate_smiles(smiles):
    # Check if the input is a valid SMILES string
    try:
        Chem.MolFromSmiles(smiles)
        return True
    except:
        return False

def process_smiles(smiles):
    # Create a molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        # Add explicit hydrogen atoms
        mol = Chem.AddHs(mol)

        # Calculate molecular weight
        mw = Descriptors.MolWt(mol)
        print("Molecular weight:", mw)

        # Display additional molecular descriptors
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        num_rings = Chem.GetSSSR(mol)  # Number of rings
        print("Number of atoms:", num_atoms)
        print("Number of bonds:", num_bonds)
        print("Number of rings:", num_rings)

        # Ask the user how they want to visualize the molecule
        visualization_option = input("How do you want to visualize the molecule? (2D/3D/Interactive): ").lower()
        
        # Draw the molecule based on the user's choice
        if visualization_option == '2d':
            Draw.MolToImage(mol).show()
        elif visualization_option == '3d':
            from rdkit.Chem import AllChem
            AllChem.EmbedMolecule(mol)
            mol_viewer = Chem.Draw.MolToMPL(mol)
            plt.show()
        elif visualization_option == 'interactive':
            import nglview
            nglview.show_rdkit(mol)
        else:
            print("Invalid visualization option. Please choose 2D, 3D, or Interactive.")
    else:
        print("Invalid SMILES string")

def main():
    # Prompt the user to enter a SMILES string
    smiles = input("Enter the SMILES string of the molecule: ")

    # Validate the SMILES input
    if validate_smiles(smiles):
        # Process the SMILES string
        process_smiles(smiles)
    else:
        print("Invalid SMILES string")

if __name__ == "__main__":
    # %matplotlib inline
    main()

#interactive molecular projection: 
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem  # Import AllChem here
import py3Dmol

def validate_smiles(smiles):
    # Check if the input is a valid SMILES string
    try:
        Chem.MolFromSmiles(smiles)
        return True
    except:
        return False

def process_smiles(smiles):
    # Create a molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        # Add explicit hydrogen atoms
        mol = Chem.AddHs(mol)

        # Calculate molecular weight
        mw = Descriptors.MolWt(mol)
        print("Molecular weight:", mw)

        # Display additional molecular descriptors
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        num_rings = Chem.GetSSSR(mol)  # Number of rings
        print("Number of atoms:", num_atoms)
        print("Number of bonds:", num_bonds)
        print("Number of rings:", num_rings)

        # Generate 3D coordinates (conformers)
        AllChem.EmbedMolecule(mol)  # Using AllChem here

        # Ask the user how they want to visualize the molecule
        visualization_option = input("How do you want to visualize the molecule? (2D/3D/Interactive): ").lower()
        
        # Draw the molecule based on the user's choice
        if visualization_option == '2d':
            Draw.MolToImage(mol).show()
        elif visualization_option == '3d':
            AllChem.EmbedMolecule(mol)  # Embed again if needed
            mol_viewer = Chem.Draw.MolToMPL(mol)
        elif visualization_option == 'interactive':
            xyz = Chem.MolToXYZBlock(mol)
            p = py3Dmol.view(width=400, height=400)
            p.addModel(xyz, 'xyz')
            p.setStyle({'stick': {}})
            p.setBackgroundColor('white')
            p.zoomTo()
            return p.show()
        else:
            print("Invalid visualization option. Please choose 2D, 3D, or Interactive.")
    else:
        print("Invalid SMILES string")

def main():
    # Prompt the user to enter a SMILES string
    smiles = input("Enter the SMILES string of the molecule: ")

    # Validate the SMILES input
    if validate_smiles(smiles):
        # Process the SMILES string
        process_smiles(smiles)
    else:
        print("Invalid SMILES string")

if __name__ == "__main__":
    main()
    
# the following code enables the 2D as well as interactive visualization at the same time of a molecule. It also displays the molecular weight, number of bonds, molecular formula and number of atoms:
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem  # Import AllChem here
import pandas as pd
import py3Dmol
from IPython.display import display

def validate_smiles(smiles):
    # Check if the input is a valid SMILES string
    try:
        Chem.MolFromSmiles(smiles)
        return True
    except:
        return False

def get_molecular_formula(mol):
    # Calculate molecular formula
    return Chem.rdMolDescriptors.CalcMolFormula(mol)

def process_smiles(smiles):
    # Create a molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        # Add explicit hydrogen atoms
        mol = Chem.AddHs(mol)

        # Calculate molecular weight
        mw = Descriptors.MolWt(mol)

        # Display additional molecular descriptors
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()

        # Embed the molecule to generate conformers
        AllChem.EmbedMolecule(mol)

        # Generate XYZ block after embedding
        mol_xyz = Chem.MolToXYZBlock(mol)
        p = py3Dmol.view(width=400, height=400)
        p.addModel(mol_xyz, 'xyz')
        p.setStyle({'stick': {}})
        p.setBackgroundColor('white')
        p.zoomTo()
        p.show()

        # Generate 2D representation and display in Jupyter Lab
        mol_img = Draw.MolToImage(mol)
        display(mol_img)

        # Get molecular formula
        molecular_formula = get_molecular_formula(mol)

        # Create DataFrame for table
        data = {
            'Property': ['Molecular Weight', 'Number of Atoms', 'Number of Bonds', 'Molecular Formula'],
            'Value': [mw, num_atoms, num_bonds, molecular_formula]
        }

        df = pd.DataFrame(data)

        # Display table
        display(df)
    else:
        print("Invalid SMILES string")

def main():
    # Prompt the user to enter a SMILES string
    smiles = input("Enter the SMILES string of the molecule: ")

    # Validate the SMILES input
    if validate_smiles(smiles):
        # Process the SMILES string
        process_smiles(smiles)
    else:
        print("Invalid SMILES string")

if __name__ == "__main__":
    main()



# FINAL CODE, displaying 2D as well as interactive 3D projection. It displays molecular characteristics (Molecular Weight, Number of Atoms, Number of Bonds, Molecular Formula, IUPAC Name, Exact Mass, Monoisotopic Mass) based on a PubChem dataset:
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import py3Dmol
from IPython.display import display

def validate_smiles(smiles):
    # Check if the input is a valid SMILES string
    try:
        Chem.MolFromSmiles(smiles)
        return True
    except:
        return False

def get_molecular_formula(mol):
    # Calculate molecular formula
    return Chem.rdMolDescriptors.CalcMolFormula(mol)

def process_smiles(smiles, df_dataset):
    # Create a molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        # Add explicit hydrogen atoms
        mol = Chem.AddHs(mol)

        # Calculate molecular weight
        mw = Descriptors.MolWt(mol)

        # Display additional molecular descriptors
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()

        # Embed the molecule to generate conformers
        AllChem.EmbedMolecule(mol)

        # Generate XYZ block after embedding
        mol_xyz = Chem.MolToXYZBlock(mol)
        p = py3Dmol.view(width=400, height=400)
        p.addModel(mol_xyz, 'xyz')
        p.setStyle({'stick': {}})
        p.setBackgroundColor('white')
        p.zoomTo()
        p.show()

        # Generate 2D representation and display in Jupyter Lab
        mol_img = Draw.MolToImage(mol)
        display(mol_img)

        # Get molecular formula
        molecular_formula = get_molecular_formula(mol)

        # Extract additional information from dataset based on canonical SMILES
        row = df_dataset.loc[df_dataset['canonicalsmiles'] == smiles].iloc[0]
        iupac_name = row['iupacname']
        exact_mass = row['exactmass']
        monoisotopic_mass = row['monoisotopicmass']

        # Create DataFrame for table
        data = {
            'Property': ['Molecular Weight', 'Number of Atoms', 'Number of Bonds', 'Molecular Formula',
                         'IUPAC Name', 'Exact Mass', 'Monoisotopic Mass'],
            'Value': [mw, num_atoms, num_bonds, molecular_formula,
                      iupac_name, exact_mass, monoisotopic_mass]
        }

        df = pd.DataFrame(data)

        # Display table
        display(df)
    else:
        print("Invalid SMILES string")

def main():
    # Prompt the user to enter a SMILES string
    smiles = input("Enter the SMILES string of the molecule: ")

    # Validate the SMILES input
    if validate_smiles(smiles):
        # Load the dataset
        dataset_path = '/path/to/your/dataset.csv'  # Replace with the actual path to your dataset
        df_dataset = pd.read_csv(dataset_path)

        # Process the SMILES string
        process_smiles(smiles, df_dataset)
    else:
        print("Invalid SMILES string")

if __name__ == "__main__":
    main()
