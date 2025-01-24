#!/usr/bin/env python3

import argparse
import os
import json
import pandas as pd
from tqdm import tqdm
from strokeDTI.predict_dti.encoder import *
from strokeDTI.predict_dti.model import *
from strokeDTI.predict_dti.data_processing import *
from strokeDTI.predict_dti.train_test_utility import *
from strokeDTI.predict_dti.samples_for_testing import *
from strokeDTI.predict_dti.params import *


def parse_arguments():
    parser = argparse.ArgumentParser(description="Compute DTI Scores")
    parser.add_argument(
        "--drug_csv",
        type=str,
        required=True,
        help="Path to the input drug CSV file (e.g., drug_list_with_smiles_first_20.csv)",
    )
    parser.add_argument(
        "--sequence_json",
        type=str,
        required=True,
        help="Path to the input sequence JSON file (e.g., sequence_dic.json)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output/",
        help="Directory to save the output DTI_output.csv (default: output/)",
    )
    parser.add_argument(
        "--total_fold", type=int, default=5, help="Number of folds to use (max: 5)"
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Load sequence dictionary
    with open(args.sequence_json, "r") as file:
        loaded_sequence_dic = json.load(file)

    model_name_list = [
        "transformer_cnn",
        "gatv2conv_cnn",
        "gineconv_cnn",
        "mpnn_cnn",
        "ResGatedGraphConv",
    ]

    # Load drug data
    drug_df = pd.read_csv(args.drug_csv)

    # Clean target sequences
    cleaned_target_dict = {
        target: remove_line(seq) for target, seq in loaded_sequence_dic.items()
    }

    # Initialize a list to store the results
    results = []

    total_iterations = (
        len(model_name_list) * args.total_fold * len(drug_df) * len(cleaned_target_dict)
    )
    with tqdm(total=total_iterations, desc="Computing DTI_scores", position=0) as pbar:
        for model_name in model_name_list:
            for fold in range(1, args.total_fold + 1):  # Assuming folds start at 1
                # Setup the model for the current fold
                test_model = setup_model(model_name, fold)

                for target, target_sequence in cleaned_target_dict.items():
                    for _, drug_row in drug_df.iterrows():
                        drug_name = drug_row["drug_names"]
                        drug_smiles = drug_row["drug_smiles"]

                        # Compute the DTI score using the test function
                        dti_score = test(test_model, drug_smiles, target_sequence)

                        # Append the result to the list
                        results.append(
                            {
                                "drug_names": drug_name,
                                "drug_smiles": drug_smiles,
                                "model": model_name,
                                "fold": fold,
                                "target": target,
                                "DTI_score": dti_score,
                            }
                        )

                        # Update the progress bar
                        pbar.update(1)

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Define output file path
    output_file = os.path.join(args.output_dir, "DTI_output.csv")

    # Save results to CSV
    results_df.to_csv(output_file, index=False)

    print(f"DTI scores have been saved to {output_file}")


if __name__ == "__main__":
    main()
