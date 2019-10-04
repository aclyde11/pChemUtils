'''
This script takes a database input and a new list. It cannonicalizes all thes miles, counts up validty and computes overlap
with the training set.
'''

import argparse

from utils import calc
from utils.files import is_valid_file_str

def main(dbase, i):
    # get valid numbers
    dbase_mols, dbase_total_mols, dbase_valid_mols = calc.pSmilesToValidMolsFromFile(dbase, threads=16,
                                                                                     return_counts=True)
    input_mols, input_total_mols, input_valid_mols = calc.pSmilesToValidMolsFromFile(i, threads=16, return_counts=True)

    print("dbase has", dbase_total_mols, "mols. Valid mols:", dbase_valid_mols)
    print("input has", input_total_mols, "mols. Valid mols:", input_valid_mols)

    # get fingprints for everything
    print("Getting fingerprints")
    dbase_fp = calc.pGetSmiles(dbase_mols, threads=4)
    input_mols = calc.pGetSmiles(input_mols, threads=4)

    # get internal consistency
    cutoff = 0.99
    print("Getting internal consistency. Cutoff:", cutoff)
    #num_repeats_inside_dbase = calc.dbaseEqualInternal(dbase_fp)
    num_repeats_inside_input = calc.dbaseEqualInternal(input_mols)
    #print("num_repeats_inside_dbase", num_repeats_inside_dbase)
    print("num_repeats_inside_input", num_repeats_inside_input)

    cutoff = 0.98
    print("Cutoff:", cutoff)
    print("Getting sim between input and database")
    counts = calc.dbaseEquality(dbase_fp, input_mols, threads=32)
    print("Total mols with >=1", counts)
    print("\n\n---------\n\n")
    print("Total Sampled: ", input_total_mols)
    print("Valid Sampled: ", input_valid_mols, float(input_valid_mols) / input_total_mols)
    print("Valid Unqiue (In Sample)", input_valid_mols - num_repeats_inside_input,
                                       float(input_valid_mols - num_repeats_inside_input) / input_total_mols)
    print("Valid Unique (W/ Trn", input_valid_mols - counts, float(input_valid_mols - counts) / input_total_mols)

# def main(dbase, i):
#     # get valid numbers
#     dbase_mols, dbase_total_mols, dbase_valid_mols = calc.pSmilesToValidMolsFromFile(dbase, threads=8,
#                                                                                      return_counts=True)
#     input_mols, input_total_mols, input_valid_mols = calc.pSmilesToValidMolsFromFile(i, threads=4, return_counts=True)
#
#     print("dbase has", dbase_total_mols, "mols. Valid mols:", dbase_valid_mols)
#     print("input has", input_total_mols, "mols. Valid mols:", input_valid_mols)
#
#     # get fingprints for everything
#     print("Getting fingerprints")
#     dbase_fp = calc.pGetFingerprints(dbase_mols, threads=4)
#     input_mols = calc.pGetFingerprints(input_mols, threads=4)
#
#     # get internal consistency
#     cutoff = 0.99
#     print("Getting internal consistency. Cutoff:", cutoff)
#     #num_repeats_inside_dbase = calc.pGetInternalSimilarity(dbase_fp, sim_cutoff=cutoff)
#     #num_repeats_inside_input = calc.pGetInternalSimilarity(input_mols, sim_cutoff=cutoff)
#     #print("num_repeats_inside_dbase", num_repeats_inside_dbase)
#     #print("num_repeats_inside_input", num_repeats_inside_input)
#
#     cutoff = 0.98
#     print("Cutoff:", cutoff)
#     print("Getting sim between input and database")
#     counts, count_total = calc.pCompareToDatabase(dbase_fp, input_mols, sim_cutoff=cutoff, threads=8)
#     print("Total mols with >=1", count_total)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-dbase', required=True, type=str, help='database list of smiles only')
    parser.add_argument('-i', required=True, type=str, help='new list of smiles only')
    args = parser.parse_args()

    main(args.dbase, args.i)
