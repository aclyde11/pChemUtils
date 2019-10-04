import multiprocessing

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

from tqdm import tqdm

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

    if return_counts:
        return iters, original_length, newlen
    return iters

def getMorganfingerprint(m):
    return AllChem.GetMorganFingerprint(m, 2)

def pGetFingerprints(mollist, threads=1):
    pool = multiprocessing.Pool(threads)
    res = pool.map(getMorganfingerprint, mollist)
    return res

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

def _pCompareDatabase(dbase, sim_cutoff=0.8):
    def _pcompare(mol):
        count = 0
        s = DataStructs.BulkTanimotoSimilarity(mol, dbase)  # +1 compare with the next to the last fp
        for m in range(len(s)):
            if s[m] >= sim_cutoff:
                count += 1

        return count

    return _pcompare

def pCompareToDatabase(dbase, mols, sim_cutoff=0.8, threads=1):
    pool = multiprocessing.Pool(threads)
    f =_pCompareDatabase(dbase, sim_cutoff)
    iters = pool.imap(f, mols, chunksize=20)

    counts = []
    count_total =0
    for count in tqdm(iters, total=len(mols)):
        counts.append(count)
        if count >= 2:
            count_total += 1

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