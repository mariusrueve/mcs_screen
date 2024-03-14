#!/usr/bin/env python
import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS
import argparse
import os
from concurrent.futures import ThreadPoolExecutor
import threading

# suppress rdkit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class MCS_Screen:
    def __init__(self, query_file, database_file, output_file, threshold):
        self.query_file = query_file
        self.database_file = database_file
        self.output_file = output_file
        self.threshold = threshold
        self.lock = threading.Lock()

        # Check if all files exist and output file does not exist
        if not os.path.isfile(self.query_file):
            raise FileNotFoundError("Query file not found")
        if not os.path.isfile(self.database_file):
            raise FileNotFoundError("Database file not found")
        if os.path.isfile(self.output_file):
            raise FileExistsError("Output file already exists")

        # check if extension is .sdf or .csv or .smi
        # if .csv or .smi assume first column is SMILES

        # get extension of query file
        self.query_ext = os.path.splitext(self.query_file)[1]
        self.database_ext = os.path.splitext(self.database_file)[1]

        print("Reading query file: ", self.query_file)
        if self.query_ext == ".sdf":
            self.query_mols = [
                mol for mol in Chem.SDMolSupplier(self.query_file) if mol is not None
            ]
        elif self.query_ext == ".csv":
            self.query_mols = [
                Chem.MolFromSmiles(line.strip().split(",")[0])
                for line in open(self.query_file)
            ]
        elif self.query_ext == ".smi":
            self.query_mols = [
                Chem.MolFromSmiles(line.strip().split()[0])
                for line in open(self.query_file)
            ]
        else:
            raise ValueError("Query file format not supported")
        
        print("Reading database file: ", self.database_file)
        if self.database_ext == ".sdf":
            self.database_mols = [
                mol for mol in Chem.SDMolSupplier(self.database_file) if mol is not None
            ]
        elif self.database_ext == ".csv":
            self.database_mols = [
                Chem.MolFromSmiles(line.strip().split(",")[0])
                for line in open(self.database_file)
            ]
        elif self.database_ext == ".smi":
            self.database_mols = [
                Chem.MolFromSmiles(line.strip().split()[0])
                for line in open(self.database_file)
            ]
        else:
            raise ValueError("Database file format not supported")

        # remove None values from the lists
        self.query_mols = [mol for mol in self.query_mols if mol is not None]
        self.database_mols = [mol for mol in self.database_mols if mol is not None]

        # self.query_mols = [
        #     mol for mol in Chem.SDMolSupplier(self.query_file) if mol is not None
        # ]
        # self.database_mols = [
        #     mol for mol in Chem.SDMolSupplier(self.database_file) if mol is not None
        # ]

    def mcs_screen(self, query_mol):
        for db_mol in self.database_mols:
            mcs = rdFMCS.FindMCS([query_mol, db_mol])
            mcs_atoms = mcs.numAtoms
            db_mol_atoms = db_mol.GetNumAtoms()
            # TODO: If there exists a small db molecule it will always results in an mcs of the size bigger than 70% of the db molecule
            # TODO: This is not a good way to filter out molecules

            # if mcs atoms or db_mol atoms is 0, skip
            if mcs_atoms == 0 or db_mol_atoms == 0:
                continue
            if mcs_atoms / db_mol_atoms >= self.threshold:
                return
        
        # write to output file SD file
        with self.lock:
            try:
                with open(self.output_file, "a") as f:
                    w = Chem.SDWriter(f)
                    w.write(query_mol)
                    w.flush()
            except Exception as e:
                print(f"Error writing to output file while processing {query_mol}: ", e)
        print("Molecule passed MCS screening: ", Chem.MolToSmiles(query_mol))

    def screen(self):
        print("Screening molecules using MCS")
        with ThreadPoolExecutor(max_workers=None) as pool:
            # pool.map(self.mcs_screen, self.query_mols)
            for _ in pool.map(self.mcs_screen, self.query_mols):
                pass
        print("MCS screening complete")


def main():
    parser = argparse.ArgumentParser(description="MCS Screening")
    parser.add_argument(
        "-q",
        "--query",
        type=str,
        required=True,
        help="Path to the query file. Accepted formats: .sdf, .csv, .smi. If .csv or .smi, first column should be SMILES."
    )
    parser.add_argument(
        "-d",
        "--database",
        type=str,
        required=True,
        help="Path to the database file. Accepted formats: .sdf, .csv, .smi. If .csv or .smi, first column should be SMILES."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        help="Path to the output file (default: mols_passed_screening.sdf)",
        default="mols_passed_screening.sdf",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        required=False,
        help="Threshold for MCS screening (default: 0.7)",
        default=0.7,
    )
    args = parser.parse_args()

    mcs_screen = MCS_Screen(args.query, args.database, args.output, args.threshold)
    mcs_screen.screen()


if __name__ == "__main__":
    main()
