from strokeDTI.predict_dti.encoder import *
from strokeDTI.predict_dti.model import get_model_from_name
from strokeDTI.predict_dti.params import DEVICE
import json


def graph_from_molecule(molecule):
    # Create a graph from a SMILE molecule that contains:
    # atom_f - Nodes and node features
    # bond_pairs - Bond pairs[atom1, atom2]
    # bond_attributes

    atom_f = []
    bond_pairs = []
    bond_attribues = []

    [atom_f.append(get_atom_features(atom)) for atom in molecule.GetAtoms()]

    for bond in molecule.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        atom1_index = atom1.GetIdx()
        atom2_index = atom2.GetIdx()

        bond_pairs.append(torch.Tensor([atom1_index, atom2_index]))
        bond_attribues.append(get_bond_features(bond))

        bond_pairs.append(torch.Tensor([atom2_index, atom1_index]))
        bond_attribues.append(get_bond_features(bond))

    atom_f = torch.stack(atom_f, 0)
    bond_pairs = torch.stack(bond_pairs, 0)
    bond_attribues = torch.stack(bond_attribues, 0)

    return atom_f, bond_pairs, bond_attribues


def get_save_path(model_name, fold_num, target):
    drug_list_save_path = "trained_model"
    save_pth = (
        drug_list_save_path
        + model_name
        + "_fold_"
        + str(fold_num)
        + "_test_target_"
        + target
        + ".csv"
    )
    return save_pth


def setup_model(modelName, foldNumber, model_output="data/trained_model/"):
    test_model = get_model_from_name(modelName).to(DEVICE)

    # load the best model checkpoint
    best_model_cp = torch.load(
        f"{model_output}best_model_" + modelName + "_fold_" + str(foldNumber) + ".pth",
        map_location=DEVICE,
        weights_only=False,
    )

    test_model.load_state_dict(best_model_cp["model_state_dict"])

    return test_model
