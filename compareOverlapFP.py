'''
This script takes a database input and a new list.

Calculates fingerprints for all compounds, and compares how many in list repeat within list, and how many are already in database.

This program streams in smiles, but has to store database and new list, so memory can be an issue.
'''

import argparse

from utils import calc
from utils.files import is_valid_file_str


def main(dbase, i):
    # get valid numbers
    dbase_mols, dbase_total_mols, dbase_valid_mols = calc.pSmilesToValidMolsFromFile(dbase, threads=8,
                                                                                     return_counts=True)
    input_mols, input_total_mols, input_valid_mols = calc.pSmilesToValidMolsFromFile(i, threads=4, return_counts=True)

    print("dbase has", dbase_total_mols, "mols. Valid mols:", dbase_valid_mols)
    print("input has", input_total_mols, "mols. Valid mols:", input_valid_mols)

    # get fingprints for everything
    print("Getting fingerprints")
    dbase_fp = calc.pGetFingerprints(dbase_mols, threads=4)
    input_mols = calc.pGetFingerprints(input_mols, threads=4)

    # get internal consistency
    print("Getting internal consistency")
    cutoff = 0.95
    num_repeats_inside_dbase = calc.pGetInternalSimilarity(dbase_fp, sim_cutoff=cutoff)
    #num_repeats_inside_input = calc.pGetInternalSimilarity(input_mols, sim_cutoff=cutoff)
    print("num_repeats_inside_dbase", num_repeats_inside_dbase)
    #print("num_repeats_inside_input", num_repeats_inside_input)

    print("Getting sim between input and database")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-dbase', required=True, type=str, help='database list of smiles only')
    parser.add_argument('-i', required=True, type=str, help='new list of smiles only')
    args = parser.parse_args()

    main(args.dbase, args.i)
