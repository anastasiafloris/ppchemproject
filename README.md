# 🧪 ppchemproject 
# The following package ✨ : 
- Validates input SMILES strings for their correctness.
- Calculates molecular weight, number of atoms, number of bonds, molecular formula.
- Generates 2D and 3D interactive representations of the molecular structure.
- Extracts additional information about the compound from a PubChem dataset (encompassing more than 100'000 molecules), giving the IUPAC name, exact mass as well as monoisotopic mass.
- The following chemical tools are used: Python, RDKit, pandas, py3Dmol, Jupyter Lab

# To install package 🖥️ :

1) Create a GitHub account at the following URL: https://github.com

2) To create a conda environment and install the required dependencies for the package:
- Run the following commands in your terminal:
- conda create -n molenv python=3.8
- conda activate molenv
- conda install -c conda-forge pandas rdkit py3Dmol ipython pytest jupyterlab

3) To activate the new environment:
- Run the following command in your terminal:
- conda activate molenv

4) To install the package: 
- Run the following command in your terminal:
- pip install git+https://github.com/anastasiafloris/ppchemproject.git

5) To navigate to the package:
- Run the following command in your terminal:
- cd ppchemproject

6) To open the package in Jupyter Lab:
- Run the following command in your terminal:
- jupyter lab

# To download the dataset 📊 :
1) Go to the following URL: https://pubchem.ncbi.nlm.nih.gov/#input_type=list&query=AUanftrxv02IZ71-Pwb0Voo9S128Ma6U1LG12M-gp9nPuZs&collection=compound&alias=PubChem%3A%20PubChem%20Compound%20TOC%3A%20Toxicity

2) Download the dataset as a CSV. file with no compression

3) Copy the dataset pathway: Navigate through Finder, right click, click on option and click on "Copy (...) as pathname".

# Run the code 🤓 :
1) Open the package module under "src/ppchemproject" and replace the dataset pathway with yours in the "moduleppchem" at line 120 where it says "Replace with the actual path to your dataset"

2) Open the Jupyter Notebook under "Notebooks" and replace the dataset pathway with yours in the Notebook titled "project_notebook" where it is indicated to "Replace with the actual path to your dataset".

3) Run the code.


