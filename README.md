# MCS_SCREEN

`mcs_screen` is a Python script that utilizes the RDKit library to screen a query file (SDF) against a database file (SDF). It calculates the maximum common substructure (MCS) between the query and the database molecules. If the size of the MCS does not exceed a user-defined threshold, the query molecule is written to the output file.

This script is particularly useful when you have a set of known inactive molecules (database) and you want to identify potentially active molecules (query) that are dissimilar to the database.

`mcs_screen` is a command-line tool that can be easily integrated into your workflow.

## Installation

```bash
git clone https://github.com/mariusrueve/mcs_screen.git
cd mcs_screen
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Usage

```bash
python mcs_screen.py -h
```

```
usage: mcs_screen.py [-h] -q QUERY -d DATABASE [-o OUTPUT] [-t THRESHOLD]

MCS Screening

options:
  -h, --help            show this help message and exit
  -q QUERY, --query QUERY
                        Path to the query file
  -d DATABASE, --database DATABASE
                        Path to the database file
  -o OUTPUT, --output OUTPUT
                        Path to the output file (default: mols_passed_screening.sdf)
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for MCS screening (default: 0.7)
```

The script accepts the following file formats as query and database files:
- SDF
- SMI (with no header and SMILES in the first column)
- CSV (with a column named "SMILES")

## Example

```bash
python mcs_screen.py --query test_data/5k.sdf --database test_data/100.sdf
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
