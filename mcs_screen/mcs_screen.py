#!/usr/bin/env python
import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS
import os
from concurrent.futures import ThreadPoolExecutor
import threading

# suppress rdkit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class mcs_screen:
    def __init__(self, query_file, database_file, output_file, threshold):
        self.query_file = query_file
        self.database_file = database_file
        self.output_file = output_file
        self.threshold = threshold
        self.lock = threading.Lock()

        # Check if all files exist and output file does not exist
        if not os.path.isfile(self.query_file):
            raise FileNotFoundError(f"Query file {self.query_file} not found")
        if not os.path.isfile(self.database_file):
            raise FileNotFoundError(f"Database file {self.database_file} not found")
        if os.path.isfile(self.output_file):
            raise FileExistsError(f"Output file {self.output_file} already exists")

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

    def screening(self, query_mol):
        for db_mol in self.database_mols:
            db_mol_atoms = db_mol.GetNumAtoms()
            if db_mol_atoms == 0:
                continue

            mcs = rdFMCS.FindMCS([query_mol, db_mol])
            mcs_atoms = mcs.numAtoms

            # if mcs atoms or db_mol atoms is 0, skip
            if mcs_atoms < 6:
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
                print(
                    f"Error writing to output file while processing \
                    {Chem.MolToSmiles(query_mol)}: ",
                    e,
                )
        print("Molecule passed MCS screening: ", Chem.MolToSmiles(query_mol))

    def start(self):
        print("Screening molecules using MCS")
        with ThreadPoolExecutor(max_workers=None) as pool:
            # pool.map(self.mcs_screen, self.query_mols)
            for _ in pool.map(self.screening, self.query_mols):
                pass
        print("MCS screening complete")
