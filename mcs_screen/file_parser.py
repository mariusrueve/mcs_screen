import os

from rdkit import Chem


def read_molecules(file):
    ext = os.path.splitext(file)[1]
    if ext == ".sdf":
        mols = [mol for mol in Chem.SDMolSupplier(file) if mol is not None]
    elif ext == ".csv":
        mols = [Chem.MolFromSmiles(line.strip().split(",")[0]) for line in open(file)]
    elif ext == ".smi":
        mols = [Chem.MolFromSmiles(line.strip().split()[0]) for line in open(file)]
    else:
        raise ValueError(f"File format {ext} not supported")

    # remove None values from the lists
    mols = [mol for mol in mols if mol is not None]
    print(f"Read {len(mols)} molecules from {file}")

    # add smiles to the molecule properties
    for mol in mols:
        mol.SetProp("SMILES", Chem.MolToSmiles(mol))

    for mol in mols:
        for atom in mol.GetAtoms():
            value = 500 * atom.GetIsAromatic() + atom.GetAtomicNum()
            atom.SetIsotope(value)

    return mols


def write_passed_mol(mol, file, writer):
    ext = os.path.splitext(file)[1]
    if ext == ".sdf":
        try:
            with open(file, "a") as _:
                writer.write(mol)
                writer.flush()
        except Exception as e:
            print(
                f"Error writing to output file while processing \
                    {mol.GetProp('SMILES')}: ",
                e,
            )
    elif ext == ".csv":
        try:
            with open(file, "a") as f:
                f.write(f"{mol.GetProp('SMILES')}\n")
        except Exception as e:
            print(
                f"Error writing to output file while processing \
                    {mol.GetProp('SMILES')}: ",
                e,
            )
    else:
        raise ValueError(f"File format {ext} not supported")


def write_failed_mol(input_mol, nonactive_mol, file):
    # if file does not exists, create it and write the header
    if not os.path.isfile(file):
        with open(file, "w") as f:
            f.write("input_smiles,nonactive_smiles\n")
    else:
        try:
            # write the input smiles and the nonactive smiles to the file
            with open(file, "a") as f:
                f.write(
                    f"{input_mol.GetProp('SMILES')},{nonactive_mol.GetProp('SMILES')}\n"
                )
        except Exception as e:
            print(
                f"Error writing to output file while processing \
                    {input_mol.GetProp('SMILES')}: ",
                e,
            )


def check_output_path(file, mode):
    if os.path.isfile(file):
        raise FileExistsError(f"Output file {file} already exists")
    else:
        if mode == "passed":
            if os.path.splitext(file)[1] == ".sdf":
                return file
            elif os.path.splitext(file)[1] == ".csv":
                return file
            else:
                raise ValueError(
                    f"Output file {file} should have .sdf or .csv extension"
                )
        elif mode == "not_passed":
            if os.path.splitext(file)[1] == ".csv":
                return file
            else:
                raise ValueError(f"Output file {file} should have .csv extension")
