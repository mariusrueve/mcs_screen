import argparse
from mcs_screen.mcs_screen import mcs_screen
import os


def main():
    parser = argparse.ArgumentParser(description="MCS Screening")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="""Path to the input file.
        Accepted formats: .sdf, .csv, .smi.
        If .csv or .smi, first column should be SMILES.""",
    )
    parser.add_argument(
        "-n",
        "--nonactives",
        type=str,
        required=True,
        help="""Path to the non-actives database file.
        Accepted formats: .sdf, .csv, .smi.
        If .csv or .smi, first column should be SMILES.""",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        help="""Path to the output file with .csv or .sdf extension
        (default: mols_passed_screening.csv)""",
        # default is input file name with _passed.sdf
        default=None,
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        required=False,
        help="""Threshold for MCS screening. Increase to filter out less molecules
        ,decrease to filter out more molecules (default: 0.7)""",
        default=0.7,
    )
    args = parser.parse_args()

    # create output file passed and not passed. If not provided, use input file name
    if args.output is None:
        passed_path = os.path.splitext(args.input)[0] + "_passed.csv"
        not_passed_path = os.path.splitext(args.input)[0] + "_not_passed.csv"
    else:
        passed_path = args.output
        not_passed_path = os.path.splitext(args.output)[0] + "_not_passed.csv"

    screen = mcs_screen(
        args.input, args.nonactives, passed_path, not_passed_path, args.threshold
    )
    screen.start()


if __name__ == "__main__":
    main()
