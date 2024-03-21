from rdkit import Chem

smiles_header_synonyms = {
    "SMILES": "SMILES",
    "smiles": "SMILES",
    "smi": "SMILES",
    "smiles": "SMILES",
    "SMI": "SMILES",
}


class InputReader:
    """
    A class for reading input files and iterating over the molecules.

    Args:
        input_file (str): The path to the input file.

    Attributes:
        input_file (str): The path to the input file.
        ext (str): The file extension of the input file.
        reader (object): The reader object used to read the input file.

    Raises:
        ValueError: If the file format is not supported.

    """

    def __init__(self, input_file):
        self.input_file = input_file
        print("Reading input file: ", self.input_file)
        self.ext = self.get_extension()
        self.reader = self.get_reader()

    def get_extension(self):
        """
        Get the file extension of the input file.

        Returns:
            str: The file extension.

        """
        return self.input_file.split(".")[-1]

    def get_reader(self):
        """
        Get the reader object based on the file extension.

        Returns:
            object: The reader object.

        Raises:
            ValueError: If the file format is not supported.

        """
        if self.ext == "sdf":
            return Chem.SDMolSupplier(self.input_file)
        elif self.ext == "csv" or self.ext == "smi":
            return open(self.input_file)
        else:
            raise ValueError(f"File format '{self.ext}' not supported")

    def add_smiles_to_properties(self, mol):
        """
        Add the SMILES string as a property to the molecule.

        Args:
            mol (object): The molecule object.

        """
        mol.SetProp("SMILES", Chem.MolToSmiles(mol))

    def mark_aromatic_atoms(self, mol):
        """
        Mark aromatic atoms in the molecule.

        Args:
            mol (object): The molecule object.

        """
        for atom in mol.GetAtoms():
            value = 500 * atom.GetIsAromatic() + atom.GetAtomicNum()
            atom.SetIsotope(value)

    def __iter__(self):
        return self

    def __next__(self):
        if self.ext == "sdf":
            mol = next(self.reader)
        elif self.ext == "csv":
            line = next(self.reader)
            smiles = line.strip().split(",")[0]
            if smiles in smiles_header_synonyms:
                return next(self)
            mol = Chem.MolFromSmiles(line.strip().split(",")[0])
        elif self.ext == "smi":
            line = next(self.reader)
            smiles = line.strip().split()[0]
            if smiles in smiles_header_synonyms:
                return next(self)
            mol = Chem.MolFromSmiles(line.strip().split()[0])
        else:
            raise ValueError(f"File format {self.ext} not supported")

        if mol is None:
            # skip None values and return the next molecule
            return next(self)

        self.add_smiles_to_properties(mol)
        self.mark_aromatic_atoms(mol)

        return mol

    def __del__(self):
        if self.ext == ".csv" or self.ext == ".smi":
            self.reader.close()

    def __str__(self):
        return f"InputReader({self.input_file})"

    def __repr__(self):
        return f"InputReader({self.input_file})"


def main():
    input_file = "test_data/45k.sdf"
    reader = InputReader(input_file)
    print(reader)
    for i in range(9):
        print(next(reader).GetProp("SMILES"))


if __name__ == "__main__":
    main()
