import argparse
from mcs_screen.mcs_screen import mcs_screen


def main():
    parser = argparse.ArgumentParser(description="MCS Screening")
    parser.add_argument(
        "-q",
        "--query",
        type=str,
        required=True,
        help="Path to the query file. Accepted formats: .sdf, .csv, .smi."
             "If .csv or .smi, first column should be SMILES.",
    )
    parser.add_argument(
        "-d",
        "--database",
        type=str,
        required=True,
        help="Path to the database file. Accepted formats: .sdf, .csv, .smi."
             "If .csv or .smi, first column should be SMILES.",
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

    screen = mcs_screen(args.query, args.database, args.output, args.threshold)
    screen.multithreading()


if __name__ == "__main__":
    main()
