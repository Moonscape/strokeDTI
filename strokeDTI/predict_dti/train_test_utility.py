from strokeDTI.predict_dti.data_processing import *
from strokeDTI.predict_dti.encoder import *
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from strokeDTI.predict_dti.params import DEVICE
import gc


def embed_drugs(smiles):

    drug_embedding = graph_from_molecule(mol_from_smiles(smiles))

    return drug_embedding


def embed_targets(seq):

    protein_embedding = protein_2_embed(trans_protein(seq))

    return protein_embedding


def create_drug_target_pairs(smiles, seq):

    drug_embedding = embed_drugs(smiles)
    protein_embedding = embed_targets(seq)

    return Data(
        x=drug_embedding[0],
        edge_index=drug_embedding[1].t().contiguous(),
        edge_attr=drug_embedding[2],
        y=[protein_embedding],
    )


def remove_line(line):
    a = line.replace("\n", "")
    a = a.replace("\t", "")
    a = a.replace(" ", "")
    return a


# def test(model, smile, sequence):
#     """
#     Function to use the model for evaluation of binding score
#     """
#     model.eval()
#     # print('Using Model for Testing')

#     data = create_drug_target_pairs(smile, sequence)

#     data_loader = DataLoader([data, data], batch_size=2)
#     # We need a DataBatch object to input into the model

#     score_holder = []

#     with torch.no_grad():
#         for i, data in enumerate(data_loader):

#             _sequence_scores = np.array(data.y)
#             _sequences = np.stack(_sequence_scores[:, 0], 0)
#             target_seq = torch.from_numpy(_sequences)

#             data = data.to(DEVICE)

#             score = model(data, target_seq)

#             score_holder.append(np.nanmean(score.cpu().numpy()))

#             # Delete variables to free memory
#             del batch, target_seq, score, score_np
#             torch.empty_cache()  # Free GPU memory

#     # Compute final mean score
#     final_score = np.nanmean(score_holder)

#     # Clean up
#     del data_loader, data, score_holder, final_score
#     torch.cuda.empty_cache()
#     gc.collect()  # Invoke garbage collector

#     return final_score


def test(model, smile, sequence):
    """
    Function to use the model for evaluation of binding score on CPU
    """
    model.eval()

    # Create data pairs
    data = create_drug_target_pairs(smile, sequence)

    # Initialize DataLoader with num_workers=0 for debugging
    data_loader = DataLoader([data, data], batch_size=2, num_workers=0)

    score_holder = []

    with torch.no_grad():
        for i, batch in enumerate(data_loader):
            # Convert batch.y to numpy array
            _sequence_scores = np.array(batch.y)
            _sequences = np.stack(_sequence_scores[:, 0], axis=0)
            target_seq = torch.from_numpy(_sequences.astype(np.float32)).to(
                DEVICE
            )  # No need to move to DEVICE since it's CPU

            # Ensure batch is on CPU
            batch = batch.to(DEVICE)

            # Forward pass
            score = model(batch, target_seq)

            # Move score to CPU and convert to numpy
            score_np = score.detach().cpu().numpy()  # Already on CPU
            mean_score = np.nanmean(score_np)
            score_holder.append(mean_score)

            # Delete intermediate variables to free memory
            del batch, target_seq, score, score_np

    # Compute final mean score
    final_score = np.nanmean(score_holder)

    # Clean up remaining variables
    del data_loader, data, score_holder
    gc.collect()  # Ensure all memory is freed

    return final_score


def test_batch(model, smiles, sequences, batch_size=32, device="cpu"):
    """
    Evaluate the model on batches of drug-sequence pairs.

    Args:
        model: The pre-trained model to evaluate.
        smiles (list): List of SMILES strings.
        sequences (list): List of protein sequences.
        batch_size (int): Number of samples per batch.
        device (str): Device to perform computation on ('cpu' or 'cuda').

    Returns:
        list: List of mean scores for each input pair.
    """
    model.eval()

    # Create a dataset and dataloader
    dataset = [
        create_drug_target_pairs(smile_single, sequence_single)
        for smile_single, sequence_single in zip(smiles, sequences)
    ]
    data_loader = DataLoader(dataset, batch_size=batch_size, num_workers=0)

    score_holder = []

    with torch.no_grad():
        for batch in data_loader:

            _sequence_scores = np.array(batch.y)  # Adjust as needed
            _sequences = np.stack(_sequence_scores[:, 0], axis=0)
            target_seq = torch.from_numpy(_sequences.astype(np.float32)).to(device)

            # Move batch data to the appropriate device
            batch_data = batch.to(device)

            # Forward pass
            scores = model(batch_data, target_seq)

            # Move scores to CPU and convert to numpy
            scores_np = scores.detach().cpu().numpy()
            score_holder.extend(scores_np)  # Assuming one score per sample

            # Clean up
            del batch_data, target_seq, scores, scores_np
            gc.collect()

    # Clean up remaining variables
    del batch_data, target_seq, scores, scores_np
    del data_loader, dataset
    gc.collect()

    return score_holder
