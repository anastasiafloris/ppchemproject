# ppchemproject
# The following package : 
Validates input SMILES strings for their correctness.
Calculates molecular weight, number of atoms, number of bonds, molecular formula.
Generates 2D and 3D interactive representations of the molecular structure.
Extracts additional information about the compound from a PubChem dataset, giving the IUPAC name, exact mass as well as monoisotopic mass. 

# To install package:
pip install git+https://github.com/anastasiafloris/ppchemproject.git
The following chemical tools are used: Python, RDKit, pandas, py3Dmol, Jupyter Lab
The dataset was downloaded at the following URL: https://pubchem.ncbi.nlm.nih.gov/. In order to insert the dataset pathway into the code, one must download the dataset, navigate through Finder, right click, click on option and click on "Copy (...) as pathname". 
To download the chemical tools one may use the following command on their terminal: pip install rdkit pandas py3Dmol
This package's code is present under "Module 1" and its final code with all the different options is from line 260 onwards. The package encompasses more than 100'000 molecules.
