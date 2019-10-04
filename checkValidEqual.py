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
    dbase_fp = calc.pGetSmiles(dbase_mols, threads=32)
    input_mols = calc.pGetSmiles(input_mols, threads=32)

    # get internal consistency
    print("Getting internal consistency.")
    #num_repeats_inside_dbase = calc.dbaseEqualInternal(dbase_fp)
    input_mols, num_repeats_inside_input = calc.dbaseEqualInternal(input_mols)
    #print("num_repeats_inside_dbase", num_repeats_inside_dbase)
    print("num_repeats_inside_input", num_repeats_inside_input)

    print("Getting sim between input and database")
    counts = calc.dbaseEquality(dbase_fp, input_mols, threads=32)
    print("Total mols with overlap", counts)

    print("Computing Scaffold Sim")
    dbase_scaffs = calc.pgetScaffsFromSmiles(dbase_fp, threads=16)
    input_scaffs = calc.pgetScaffsFromSmiles(input_mols, threads=16)

    dbase_scaffs, _ = calc.dbaseEqualInternal(dbase_scaffs)
    input_scaffs, num_repeated_scaffs = calc.dbaseEqualInternal(input_scaffs)


    scaff_counts = calc.dbaseEquality(dbase_scaffs, input_scaffs, threads=32)


    print("\n\n---------\n\n")
    print("Total Sampled: ", input_total_mols)
    print("Valid Sampled: ", input_valid_mols, float(input_valid_mols) / input_total_mols)
    print("Valid Unqiue (In Sample)", len(input_mols),
                                       float(len(input_mols)) / input_total_mols)
    print("Valid Unique (W/ Trn", len(input_mols) - counts, float(len(input_mols) - counts) / input_total_mols)

    print("Unqiue Scaffs (In Sample)", len(input_scaffs), float(len(input_scaffs)) / float(len(dbase_scaffs)))
    print("Unique Scaffs (Inn Training", len(input_scaffs) - counts)

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
