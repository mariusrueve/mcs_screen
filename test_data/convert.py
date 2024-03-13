"""Convert all .sdf files in the current directory to .csv and .smi files, where the first column is SMILES."""
import os
from rdkit import Chem
import pandas as pd

# get all .sdf files in the current directory
sdf_files = [f for f in os.listdir() if f.endswith(".sdf")]

# convert all .sdf files to .csv and .smi
for sdf_file in sdf_files:
    # read .sdf file
    mols = [mol for mol in Chem.SDMolSupplier(sdf_file) if mol is not None]
    # convert to .csv
    df = pd.DataFrame([Chem.MolToSmiles(mol) for mol in mols], columns=["SMILES"])
    df.to_csv(sdf_file.replace(".sdf", ".csv"), index=False)
    # convert to .smi with \t as delimiter and a molecule name
    mol_names = [f"mol{i}" for i in range(len(mols))]
    with open(sdf_file.replace(".sdf", ".smi"), "w") as f:
        for mol, mol_name in zip(mols, mol_names):
            f.write(f"{Chem.MolToSmiles(mol)}\t{mol_name}\n")
    print(f"Converted {sdf_file} to .csv and .smi")
    