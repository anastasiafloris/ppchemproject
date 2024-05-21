
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
