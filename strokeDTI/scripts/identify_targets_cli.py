import argparse
import os
from strokeDTI.target_identification.pathway_merge import (
    combine_dataframes,
    combine_kegg_rna_seq,
)
from strokeDTI.target_identification.plot_merged_graph import (
    make_polar_plot,
    plot_graph,
)
import plotly.io as pio
from pathlib import Path


def read_kegg(file_path):
    with open(file_path, "r") as file:
        content = file.read()
    number_list = [num.strip() for num in content.split(",")]
    return number_list


def main():
    parser = argparse.ArgumentParser(
        description="Identify targets based on KEGG pathways and RNA-Seq data."
    )

    # Define command-line arguments
    parser.add_argument(
        "-s",
        "--species",
        required=True,
        choices=["mouse", "rat", "human"],  # Add other supported species
        help="Species name (e.g., mouse, rat, human)",
    )
    parser.add_argument(
        "-rna_seq",
        required=True,
        help="Path to RNA-Seq data file (e.g., ../data/galaxy_mouse.tabular)",
    )
    parser.add_argument(
        "-kegg",
        required=True,
        help="Path to KEGG pathways list file (e.g., kegg_list.txt)",
    )
    parser.add_argument(
        "-out",
        "--output",
        required=True,
        help="Output root directory path (e.g., ../output/)",
    )

    args = parser.parse_args()

    # Validate paths
    if not os.path.exists(args.rna_seq):
        parser.error(f"RNA-Seq file not found: {args.rna_seq}")
    if not os.path.exists(args.kegg):
        parser.error(f"KEGG list file not found: {args.kegg}")
    if not os.path.isdir(args.output):
        os.makedirs(args.output, exist_ok=True)

    # Read KEGG pathways from the provided file
    kegg_path_list = read_kegg(args.kegg)

    root_path = os.path.join(args.output, "")
    species = args.species
    rna_seq_path = args.rna_seq

    # Execute functions
    combine_dataframes(
        root_path=root_path, pathway_list=kegg_path_list, species=species, save=True
    )

    combined_df = combine_kegg_rna_seq(
        root_path=root_path,
        rna_seq_path=rna_seq_path,
        right_on="Gene name",
        species=species,
    )

    polor_plot1 = make_polar_plot(combined_df)

    image_path = root_path + "graph1.png"
    pio.write_image(polor_plot1, image_path, scale=2)

    plot_graph(
        root_path=root_path,
        kegg_path_list=kegg_path_list,
        species=species,
        sequence_path=rna_seq_path,
    )

    print(f"Analysis complete. Results saved in {root_path}")


if __name__ == "__main__":
    main()
