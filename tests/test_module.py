import pytest
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# Import the functions from the module you want to test
from ppchemproject.moduleppchem import validate_smiles, get_molecular_formula, process_smiles

def test_validate_smiles():
    # Test valid SMILES strings
    assert validate_smiles('CCO') == True
    assert validate_smiles('C1=CC=CC=C1') == True
    
    # Test invalid SMILES strings
    assert validate_smiles('helloooo') == False
    assert validate_smiles('invalid_smiles') == False

def test_get_molecular_formula():
    mol = Chem.MolFromSmiles('CCO')
    assert get_molecular_formula(mol) == 'C2H6O'
    
    mol = Chem.MolFromSmiles('C1=CC=CC=C1')
    assert get_molecular_formula(mol) == 'C6H6'

def test_process_smiles(monkeypatch):
    # Create a mock dataset
    data = {
        'canonicalsmiles': ['CCO', 'C1=CC=CC=C1'],
        'iupacname': ['ethanol', 'benzene'],
        'exactmass': [46.041864812, 78.046950192],
        'monoisotopicmass': [46.041864812, 78.046950192]
    }
    df_dataset = pd.DataFrame(data)

    # Mock the display function to prevent actual display in tests
    def mock_display(*args, **kwargs):
        pass

    monkeypatch.setattr("IPython.display.display", mock_display)
    
    # Mock the input function to provide a fixed SMILES string
    monkeypatch.setattr('builtins.input', lambda _: 'CCO')

    # Mock the py3Dmol.view function to prevent actual rendering
    class MockView:
        def __init__(self, *args, **kwargs):
            pass

        def addModel(self, *args, **kwargs):
            pass

        def setStyle(self, *args, **kwargs):
            pass

        def setBackgroundColor(self, *args, **kwargs):
            pass

        def zoomTo(self, *args, **kwargs):
            pass

        def show(self, *args, **kwargs):
            pass

    monkeypatch.setattr('py3Dmol.view', MockView)

    # Test with a valid SMILES string
    process_smiles('CCO', df_dataset)
    # Verify expected outputs, values would have to be manually checked from function calls
    
    # Test with an invalid SMILES string
    with pytest.raises(ValueError):
        process_smiles('invalid_smiles', df_dataset)

if __name__ == "__main__":
    pytest.main()
