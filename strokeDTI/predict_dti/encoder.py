import numpy as np
import torch
from rdkit import Chem
from sklearn.preprocessing import OneHotEncoder

# Basic Configurations
ELEM_LIST = [
    "C",
    "N",
    "O",
    "S",
    "F",
    "Si",
    "P",
    "Cl",
    "Br",
    "Mg",
    "Na",
    "Ca",
    "Fe",
    "Al",
    "I",
    "B",
    "K",
    "Se",
    "Zn",
    "H",
    "Cu",
    "Mn",
    "unknown",
]

RESIDUE_LIST = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLU",
    "GLN",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]


# Utility Functions
def mol_from_smiles(smiles):
    # Get the molecule from SMILE strings
    molecule = Chem.MolFromSmiles(smiles)

    # Kekulization
    Chem.Kekulize(molecule)

    if molecule == None:
        raise Exception("No Molecule")
    return molecule


def one_hot_encode(input, categories):

    holder = np.zeros(len(categories))

    if input in categories:

        holder[[i for i, x in enumerate(categories) if x == input]] = 1

    else:
        holder[-1] = 1.0

    return list(holder)


def numerical_encode(input, categories):
    if input not in categories:
        return [0]
    else:
        return [input]


def get_atom_features(atom):

    return torch.Tensor(
        one_hot_encode(atom.GetSymbol(), ELEM_LIST)
        + numerical_encode(atom.GetDegree(), [0, 1, 2, 3, 4, 5, 6])
        + numerical_encode(atom.GetFormalCharge(), [-2, -1, 0, 1, 2])
        + numerical_encode(int(atom.GetChiralTag()), [0, 1, 2, 3])
        + numerical_encode(atom.GetTotalNumHs(), [0, 1, 2, 3, 4])
        + one_hot_encode(
            atom.GetHybridization().name.lower(), ["s", "sp", "sp2", "sp3"]
        )
        + numerical_encode(atom.GetTotalValence(), [0, 1, 2, 3, 4, 5, 6])
        + [atom.GetIsAromatic()]
    )


def get_bond_features(bond):
    stereo = int(bond.GetStereo())

    return torch.Tensor(
        one_hot_encode(
            bond.GetBondType().name.lower(), ["single", "double", "triple", "aromatic"]
        )
        + numerical_encode(stereo, [0, 1, 2, 3, 4, 5])
        + [bond.GetIsConjugated()]
        + [bond.IsInRing()]
    )


# Method to transform ammino acid sequences to encodings

# NEED TO CHANGE THIS#


amino_char = [
    "?",
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "M",
    "L",
    "O",
    "N",
    "Q",
    "P",
    "S",
    "R",
    "U",
    "T",
    "W",
    "V",
    "Y",
    "X",
    "Z",
]

MAX_SEQ_PROTEIN = 1000


def protein_to_features(sequence):
    temp = list(sequence.upper())
    temp = [i if i in amino_char else "?" for i in temp]
    if len(temp) < MAX_SEQ_PROTEIN:
        temp = temp + ["?"] * (MAX_SEQ_PROTEIN - len(temp))
    else:
        temp = temp[:MAX_SEQ_PROTEIN]
    return temp


def trans_protein(x):
    temp = list(x.upper())
    temp = [i if i in amino_char else "?" for i in temp]
    if len(temp) < MAX_SEQ_PROTEIN:
        temp = temp + ["?"] * (MAX_SEQ_PROTEIN - len(temp))
    else:
        temp = temp[:MAX_SEQ_PROTEIN]
    return temp


enc_protein = OneHotEncoder().fit(np.array(amino_char).reshape(-1, 1))


def protein_2_embed(x):
    return enc_protein.transform(np.array(x).reshape(-1, 1)).toarray().T
