import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import py3Dmol
from IPython.display import display

def validate_smiles(smiles: str) -> bool:
    """
    Check if the input is a valid SMILES string.

    Parameters
    ----------
    smiles : str
        A text string representing a SMILES notation.

    Returns
    -------
    bool
        True if the input is a valid SMILES string, False otherwise.
    """
    try:
        Chem.MolFromSmiles(smiles)
        return True
    except:
        return False

def get_molecular_formula(mol):
    """
    Calculate molecular formula.

    Parameters
    ----------
    mol : RDKit Mol object
        RDKit Mol object representing the molecule.

    Returns
    -------
    str
        Molecular formula of the molecule.
    """
    return Chem.rdMolDescriptors.CalcMolFormula(mol)

def process_smiles(smiles, df_dataset):
    """
    Process the SMILES string and display molecular information.

    Parameters
    ----------
    smiles : str
        A text string representing a SMILES notation.
    df_dataset : pandas DataFrame
        DataFrame containing additional information about molecules.

    Returns
    -------
    None
    """
    # Create a molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        # Add explicit hydrogen atoms
        mol = Chem.AddHs(mol)

        # Calculate molecular weight
        mw = Descriptors.MolWt(mol)

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

        # Get molecular formula
        molecular_formula = get_molecular_formula(mol)

        # Extract additional information from dataset based on canonical SMILES
        row = df_dataset.loc[df_dataset['canonicalsmiles'] == smiles].iloc[0]
        iupac_name = row['iupacname']
        exact_mass = row['exactmass']
        monoisotopic_mass = row['monoisotopicmass']

        # Display table
        display({
            'Property': ['Molecular Weight', 'Molecular Formula', 'IUPAC Name', 'Exact Mass', 'Monoisotopic Mass'],
            'Value': [mw, molecular_formula, iupac_name, exact_mass, monoisotopic_mass]
        })
    else:
        print("Invalid SMILES string")

def main():
    """
    Main function to process SMILES input and display molecular information.
    """
    smiles = input("Enter the SMILES string of the molecule: ")
    # Load the dataset
    dataset_path = '/path/to/dataset.csv'
    df_dataset = pd.read_csv(dataset_path)
    # Process the SMILES string
    if validate_smiles(smiles):
        process_smiles(smiles, df_dataset)
    else:
        print("Invalid SMILES string")

if __name__ == "__main__":
    main()
