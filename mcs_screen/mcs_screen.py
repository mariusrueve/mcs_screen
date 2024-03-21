#!/usr/bin/env python
import threading
from concurrent.futures import ThreadPoolExecutor

import rdkit
from mcs_screen.file_parser import (
    read_molecules,
    write_passed_mol,
    write_failed_mol,
    check_output_path,
)
from mcs_screen.input_reader import InputReader

from rdkit import Chem
from rdkit.Chem import rdFMCS

# suppress rdkit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class mcs_screen:
    def __init__(
        self, input_file, nonactives_file, passed_file, not_passed_file, threshold
    ):
        self.input_file = input_file
        self.nonactives_file = nonactives_file
        self.passed_file = passed_file
        self.not_passed_file = not_passed_file

        # check if output files exists
        check_output_path(self.passed_file, "passed")
        check_output_path(self.not_passed_file, "not passed")

        self.passed_file_writer = Chem.SDWriter(self.passed_file)

        self.threshold = threshold
        self.lock = threading.Lock()

        self.input_mols_reader = InputReader(self.input_file)
        self.nonactive_mols = []

        print("Reading database file: ", self.nonactives_file)
        self.nonactive_mols = read_molecules(self.nonactives_file)

    def screening(self):
        with self.lock:
            input_mol = next(self.input_mols_reader)
            print("Screening molecule: ", input_mol.GetProp("SMILES"))
        for nonactive in self.nonactive_mols:
            input_atoms = input_mol.GetNumAtoms()
            nonactive_atoms = nonactive.GetNumAtoms()
            if nonactive_atoms == 0 or input_atoms == 0:
                continue

            mcs = rdFMCS.FindMCS(
                [input_mol, nonactive],
                ringMatchesRingOnly=True,
                completeRingsOnly=False,
                atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
                bondCompare=rdFMCS.BondCompare.CompareAny,
                timeout=60,
            )
            mcs_atoms = mcs.numAtoms

            # if mcs atoms or db_mol atoms is 0, skip
            if mcs_atoms < 6:
                continue
            if (
                nonactive_atoms > input_atoms
                and mcs_atoms / nonactive_atoms >= self.threshold
            ):
                with self.lock:
                    write_failed_mol(
                        input_mol,
                        nonactive,
                        self.not_passed_file,
                    )
                return
            elif (
                input_atoms > nonactive_atoms
                and mcs_atoms / input_atoms >= self.threshold
            ):
                with self.lock:
                    write_failed_mol(
                        input_mol,
                        nonactive,
                        self.not_passed_file,
                    )
                return

        # write to output file SD file
        with self.lock:
            write_passed_mol(input_mol, self.passed_file, self.passed_file_writer)
        print("Molecule passed MCS screening: ", input_mol.GetProp("SMILES"))

    def start(self):
        print("Screening molecules using MCS")
        with ThreadPoolExecutor(max_workers=None) as pool:
            # start screening with multiple threads
            while True:
                try:
                    pool.submit(self.screening)
                except StopIteration:
                    break
        print("MCS screening complete")
