import multiprocessing

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem


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
    iter = pool.imap(smiles, chunksize=1000)
    iter = list(filter(lambda x: x is not None, iter))
    newlen = len(iter)

    if return_counts:
        return iter, original_length, newlen
    return iter


def pGetFingerprints(mollist, threads=1):
    func = lambda m: AllChem.GetMorganFingerprint(m, 2)
    pool = multiprocessing.Pool(threads)
    res = pool.map(func, mollist)
    return res


def pGetInternalSimilarity(fps, sim_cutoff=0.8):
    # the list for the dataframe
    count = 0
    # compare all fp pairwise without duplicates
    for n in range(len(fps) - 1):  # -1 so the last fp will not be used
        s = DataStructs.BulkTanimotoSimilarity(fps[n], fps[n + 1:])  # +1 compare with the next to the last fp
        # collect the SMILES and values
        for m in range(len(s)):
            if s[m] >= sim_cutoff:
                count += 1
    return count
