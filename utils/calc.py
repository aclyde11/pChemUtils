import multiprocessing

import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import AllChem
from tqdm import tqdm

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

def noneCheckedSmileToMol(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        return mol
    except KeyboardInterrupt:
        print("Caught. Exiting")
        exit()
    except:
        return None


def FromTxtFileToSmilesList(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
        lines = list(filter(lambda x: len(x) != 0, map(lambda x: x.strip(), lines)))
    return lines


def pSmilesToValidMolsFromFile(fname, threads=1, return_counts=True):
    pool = multiprocessing.Pool(threads)

    smiles = FromTxtFileToSmilesList(fname)
    original_length = len(smiles)
    iters = pool.imap(noneCheckedSmileToMol, smiles, chunksize=1000)
    iters = list(filter(lambda x: x is not None, iters))
    newlen = len(iters)
    pool.close()

    if return_counts:
        return iters, original_length, newlen
    return iters


def getMorganfingerprint(m):
    return AllChem.GetMorganFingerprint(m, 2)


def pGetFingerprints(mollist, threads=1):
    pool = multiprocessing.Pool(threads)
    res = pool.map(getMorganfingerprint, mollist)
    pool.close()
    return res


def getSmiles(m):
    return Chem.MolToSmiles(m)


def pGetSmiles(mollist, threads=1):
    pool = multiprocessing.Pool(threads)
    res = pool.map(getSmiles, mollist)

    pool.close()
    return res


def check_in(ins):
    count = 0
    mols, dbase = ins
    for mol in mols:
        if mol in dbase:
            count += 1

    return count


def dbaseEquality(dbase, mols, threads=16):
    count = 0
    count = len(set(dbase).intersection(set(mols)))
    return count

    # mols =
    #
    # mols = np.array_split(np.array(mols), threads)
    # inputs = map(lambda x: (list(mols[x]), dbase), range(threads))
    # pool = multiprocessing.Pool(threads)
    # res = pool.map(check_in, inputs)
    # pool.close()
    #
    # sum = 0
    # for i in res:
    #     sum += i
    # return sum


def dbaseEqualInternal(dbase):
    x = set(dbase)
    return list(x), len(x)


def pGetInternalSimilarity(fps, sim_cutoff=0.8):
    # the list for the dataframe
    count = 0
    # compare all fp pairwise without duplicates
    for n in tqdm(range(len(fps) - 1)):  # -1 so the last fp will not be used
        s = DataStructs.BulkTanimotoSimilarity(fps[n], fps[n + 1:])  # +1 compare with the next to the last fp
        # collect the SMILES and values
        for m in range(len(s)):
            if s[m] >= sim_cutoff:
                count += 1
    return count


def _pCompareDatabase(ins, sim_cutoff=0.99):
    count = 0
    mol, dbase = ins
    s = DataStructs.BulkTanimotoSimilarity(mol, dbase)  # +1 compare with the next to the last fp
    for m in range(len(s)):
        if s[m] >= sim_cutoff:
            count += 1

    return count


def pCompareToDatabase(dbase, mols, sim_cutoff=0.8, threads=1):
    pool = multiprocessing.Pool(threads)

    ins = map(lambda x: (x, dbase), mols)
    iters = pool.imap(_pCompareDatabase, ins, chunksize=20)

    counts = []
    count_total = 0
    for count in tqdm(iters, total=len(mols)):
        counts.append(count)
        if count >= 2:
            count_total += 1
    pool.close()

    return counts, count_total


def compareToDatabase(dbase, mols, sim_cutoff=0.8):
    count_total = 0
    counts = []
    for mol in tqdm(mols):
        count = 0
        s = DataStructs.BulkTanimotoSimilarity(mol, dbase)  # +1 compare with the next to the last fp
        for m in range(len(s)):
            if s[m] >= sim_cutoff:
                count += 1
        counts.append(count)
        if count >= 1:
            count_total += 1
    return counts, count_total
