# Standard library packages
import io
import os
import os.path
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
import requests
import pandas as pd
from tqdm import tqdm
import re
from keggx.keggx import KEGG
import numpy as np
import networkx as nx
import logging
from pathlib import Path


def get_kegg_term_from_ko(term, species):
    x = REST.kegg_get(term).read()

    if species == "human":
        prefx = "hsa"
    elif species == "mouse":
        prefx = "mmu"
    elif species == "rat":
        prefx = "rno"

    x = x.split(prefx.upper() + ":")[1]
    x = x.split("(")[0]
    x = x.strip()
    return prefx + ":" + x


def get_uniprot_id(code_id):
    try:
        url = f"https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cgene_names%2Corganism_name&format=tsv&query=%28{code_id}%29"
        all_fastas = requests.get(url).text

        uni_id = all_fastas.split()[6]
        # print(f'Found uni_prot id {uni_id} for ensemble id {code_id}')

        if len(uni_id) == 0:
            raise ValueError("Nothing found for this target")

        else:
            return uni_id

    except:
        try:
            url = f"https://rest.uniprot.org/uniprotkb/stream?fields=accession&format=tsv&query=%28{code_id}%29"
            all_fastas = requests.get(url).text

            uni_id = all_fastas.split()[6]
            # print(f'Found uni_prot id {uni_id} for ensemble id {code_id}')

            if len(uni_id) == 0:
                raise ValueError("Nothing found for this target")

            else:
                return uni_id

        except:
            print(f"{code_id} is an invalid code")


def add_uniprot_to_df(dataframe):
    uni_prot_holder = []
    for idx in tqdm(dataframe.index):
        ens = dataframe.at[idx, "Ens"]
        uni_prot_holder.append(get_uniprot_id(ens))

    dataframe["uniprot"] = uni_prot_holder

    return dataframe


def to_df(result):
    return pd.read_table(io.StringIO(result), header=None)


def remove_line(line):
    a = line.replace("\n", "")
    a = a.replace("\t", "")
    return a


def find_sequences(text):
    x = text.split("\n")[1:]
    return "".join(x)


def get_symbol(text):
    try:
        _x = text.split("SYMBOL")[1]
        _x = _x.split("\n")[0]
        return _x.split(",")[0].strip()
    except:
        return "No SYMBOL ID avalible"


def find_uniprot(text):
    if "UniProt" in text:
        _y = text.split("UniProt")[1]
        _y = _y.split("n")[0]
        return re.sub(r"\W+", "", _y)
    else:
        return "No Uniprot ID"


def get_pathway_name(mapID):
    _pname = REST.kegg_get(mapID).read()
    _pname = _pname.split("NAME")[1]
    _pname = _pname.split("-")[0]
    _pname = _pname.strip()
    return _pname.replace(" ", "_")


def list_from_kegg_ko_entry(entry, species):
    if species == "human":
        re_pattern = r"ENSG0\w+"
    elif species == "mouse":
        re_pattern = r"ENSMUSG0\w+"
    elif species == "rat":
        re_pattern = r"ENSRNOG0\w+"
    else:
        raise AssertionError("Unknown mapping species")

    def find_ens(text):
        x = re.findall(re_pattern, text)
        if len(x) == 1:
            return x[0]
        elif len(x) == 2:
            return x[0]
        elif len(x) == 0:
            return "No Ens Found"
        else:
            raise AssertionError(
                f"Found {len(x)} ensembl number, More than one Ensembl number!"
            )

    ens_list = []
    sequence_list = []
    name_list = []
    uniprot_list = []

    term = get_kegg_term_from_ko(entry, species)

    _x = REST.kegg_get(term).read()

    ens_list.append(find_ens(_x))
    name_list.append(get_symbol(_x))
    uniprot_list.append(find_uniprot(_x))

    result = REST.kegg_get(term, "aaseq").read()
    sequence_list.append(find_sequences(result))

    # return ens_list, sequence_list, name_list, uniprot_list
    return find_ens(_x), find_sequences(result), get_symbol(_x)


def get_ko_terms_from_df(custom_df, global_species):
    ens_list = []
    sequence_list = []
    name_list = []

    for idx in custom_df.index:
        terms = custom_df.at[idx, "ko_terms"]

        ens, seq, sym = list_from_kegg_ko_entry(terms, global_species)

        ens_list.append(ens)
        sequence_list.append(seq)
        name_list.append(sym)

    custom_df["Ens"] = ens_list
    custom_df["sequences"] = sequence_list
    custom_df["name"] = name_list

    return custom_df


def df_from_kegg_map(mapID):

    if mapID[:3] == "hsa":
        mapping_id = "hsa"
        re_pattern = r"ENSG0\w+"
    elif mapID[:3] == "mmu":
        mapping_id = "mmu"
        re_pattern = r"ENSMUSG0\w+"
    elif mapID[:3] == "rno":
        mapping_id = "rno"
        re_pattern = r"ENSRNOG0\w+"
    else:
        raise AssertionError("Unknown mapping species")

    def find_ens(text):
        x = re.findall(re_pattern, text)
        if len(x) == 1:
            return x[0]
        elif len(x) == 2:
            return x[0]
        elif len(x) == 0:
            return "No Ens Found"
        else:
            raise AssertionError(
                f"Found {len(x)} ensembl number, More than one Ensembl number!"
            )

    result = REST.kegg_link(mapping_id, mapID).read()
    df = to_df(result)

    ens_list = []
    sequence_list = []
    name_list = []
    uniprot_list = []

    for index, value in tqdm(enumerate(df.iloc[:, 1])):
        result = REST.kegg_get(value).read()
        ens_list.append(find_ens(result))
        name_list.append(get_symbol(result))
        uniprot_list.append(find_uniprot(result))

        result = REST.kegg_get(value, "aaseq").read()
        sequence_list.append(find_sequences(result))

    df["Ens"] = ens_list
    df["sequences"] = sequence_list
    df["name"] = name_list

    return df


def save_dataframe_from_path_list(root_path, path_list, species):
    """
    Saves dataframes for each pathway in the path_list to CSV files.
    Skips pathways if the corresponding CSV file already exists.

    Parameters:
        root_path (str or Path): The directory where CSV files will be saved.
        path_list (list of str): List of pathway identifiers.
        species (str): Species name ('human', 'mouse', or 'rat').

    Raises:
        ValueError: If the species is not one of the expected values.
    """
    # Set up logging
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    logging.info(
        "Retrieving pathway data from KEGG and saving to CSV files. Takes around 4min per pathway."
    )

    # Mapping species to their prefixes
    species_prefixes = {"human": "hsa", "mouse": "mmu", "rat": "rno"}

    try:
        prefix = species_prefixes[species.lower()]
    except KeyError:
        raise ValueError(
            f"Unsupported species '{species}'. Supported species are: {', '.join(species_prefixes.keys())}."
        )

    root = Path(root_path)

    for identifier in path_list:
        prefixed_id = f"{prefix}{identifier}"
        pathway_name = get_pathway_name(prefixed_id)
        filename = f"{pathway_name}_{species.lower()}.csv"
        file_path = root / filename

        if file_path.exists():
            logging.info(f"File '{filename}' already exists. Skipping.")
            continue

        try:
            df = df_from_kegg_map(prefixed_id)
            logging.info(f"Saving pathway '{pathway_name}' to '{filename}'.")
            df.to_csv(file_path)
        except Exception as e:
            logging.error(f"Failed to save pathway '{pathway_name}': {e}")


def combine_dataframes(root_path, pathway_list, species, save=False):
    """
    Combines multiple pathway-specific CSV files into a single dataframe.
    Optionally saves individual and combined dataframes to CSV files.

    Parameters:
        root_path (str or Path): The directory where CSV files are located and/or will be saved.
        pathway_list (list of str): List of pathway identifiers.
        species (str): Species name ('human', 'mouse', or 'rat').
        save (bool): If True, saves individual dataframes and the combined dataframe to CSV files.

    Returns:
        pd.DataFrame: The combined dataframe containing all pathways.

    Raises:
        ValueError: If the species is not supported.
        Exception: For any unexpected errors during file operations.
    """
    # Set up logging
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    # Mapping species to their prefixes
    species_prefixes = {"human": "hsa", "mouse": "mmu", "rat": "rno"}

    species_key = species.lower()
    prefix = species_prefixes.get(species_key)
    if not prefix:
        raise ValueError(
            f"Unsupported species '{species}'. Supported species are: {', '.join(species_prefixes.keys())}."
        )

    # Ensure root_path is a Path object and exists
    root = Path(root_path)
    if not root.exists():
        if save:
            try:
                root.mkdir(parents=True, exist_ok=True)
                logging.info(f"Created root directory: {root}")
            except Exception as e:
                logging.error(f"Failed to create root directory '{root}': {e}")
                raise
        else:
            logging.warning(
                f"Root directory '{root}' does not exist. Proceeding without saving."
            )

    dataframe_list = []

    for identifier in pathway_list:
        prefixed_id = f"{prefix}{identifier}"
        pathway_name = get_pathway_name(prefixed_id)
        csv_filename = f"{pathway_name}_{species_key}.csv"
        csv_path = root / csv_filename

        try:
            if csv_path.exists():
                logging.info(f"Reading existing CSV file: '{csv_filename}'.")
                df = pd.read_csv(csv_path, index_col=0)
            else:

                logging.info(
                    f"CSV file '{csv_filename}' not found. Generating dataframe from KEGG map."
                )

                logging.info(
                    "Retrieving pathway data from KEGG and saving to CSV files. Takes around 7 min per pathway."
                )

                df = df_from_kegg_map(prefixed_id)
                if save:
                    df.to_csv(csv_path)
                    logging.info(f"Saved dataframe to '{csv_filename}'.")

            # Add pathway information
            df["pathway"] = pathway_name
            dataframe_list.append(df)

        except Exception as e:
            logging.error(f"Error processing pathway '{pathway_name}': {e}")
            continue  # Skip to the next pathway

    if not dataframe_list:
        logging.warning("No dataframes were combined. Returning an empty dataframe.")
        return pd.DataFrame()

    try:
        # Combine all dataframes
        df_combined = pd.concat(dataframe_list, ignore_index=True)

        # Rename columns if they exist
        rename_dict = {}
        if len(df_combined.columns) > 0:
            rename_dict[df_combined.columns[0]] = "pathway/gene_id"
        if len(df_combined.columns) > 1:
            rename_dict[df_combined.columns[1]] = "ko_terms"

        if rename_dict:
            df_combined.rename(columns=rename_dict, inplace=True)
            logging.info("Renamed dataframe columns for clarity.")

        # Save the combined dataframe
        if save:
            combined_filename = f"combined_pathway_{species_key}.csv"
            combined_path = root / combined_filename
            df_combined.to_csv(combined_path, index=False)
            logging.info(f"Saved combined dataframe to '{combined_filename}'.")

    except Exception as e:
        logging.error(f"Failed to combine dataframes: {e}")
        raise


def combine_kegg_rna_seq(
    root_path, rna_seq_path, species, left_on="name", right_on="GeneID"
):

    kegg_path = root_path + "combined_pathway" + "_" + species + ".csv"

    kegg_df = pd.read_csv(kegg_path, index_col=0)
    rna_seq_df = pd.read_csv(rna_seq_path, delimiter="\t")

    merged_df = pd.merge(
        kegg_df, rna_seq_df, how="inner", left_on=left_on, right_on=right_on
    )
    merged_df["fold_change"] = merged_df["log2(FC)"].apply(lambda x: 2**x)

    return merged_df


def get_kegg_list(plist, species):
    """
    Generate a list of KEGG identifiers based on the input list and species.

    Args:
        plist (List[str]): A list of input identifiers.
        species (str): The species for which to generate KEGG identifiers.

    Returns:
        List[str]: A list of KEGG identifiers generated based on the input list and species.
    """
    species_prefix = {"human": "hsa", "mouse": "mmu", "rat": "rno"}
    prefix = species_prefix.get(species, "")

    input_list = [prefix + i for i in plist]
    return input_list


def get_whole_graph(pathway_id):
    pathway = KEGG(pathway_id=pathway_id)
    pathway_edges = pathway._get_directed_edge_attributes_as_dataframe(
        pathway.edge_attributes_df
    )

    pathway_entries = pathway.entry_attributes_df
    pathway_entries = pathway_entries.replace("", np.nan).dropna(subset=["name"])

    g = nx.from_pandas_edgelist(
        pathway_edges, "source", "target", edge_attr=True, create_using=nx.DiGraph()
    )

    nx.set_node_attributes(g, pathway_entries.to_dict("index"))

    nx.relabel_nodes(g, {i: g.nodes[i]["name"] for i in g.nodes}, copy=False)

    return g


def merge_pathway_with_seq(entry_df, seq_df):
    return pd.merge(
        entry_df[entry_df.type == "gene"].reset_index(),
        seq_df,
        how="inner",
        left_on="name",
        right_on="Gene name",
    )


def select_only_connected(merge_dataframe, edges_dataframe):

    _node_list = merge_dataframe["id"]

    edges_selected = edges_dataframe[
        (edges_dataframe["source"].isin(_node_list))
        & (edges_dataframe["target"].isin(_node_list))
    ]

    connect_nodes = list(edges_selected.source) + list(edges_selected.target)
    connect_nodes = list(set(connect_nodes))

    return (
        merge_dataframe[merge_dataframe["id"].isin(connect_nodes)].set_index("id"),
        edges_selected,
    )


def relabel_nodes(graph):
    nx.relabel_nodes(
        graph, {i: graph.nodes[i]["name_x"] for i in graph.nodes}, copy=False
    )


def get_graph_from_df(pathway_id, seq_df):
    pathway = KEGG(pathway_id=pathway_id)
    pathway_edges = pathway._get_directed_edge_attributes_as_dataframe(
        pathway.edge_attributes_df
    )

    pathway_entries = pathway.entry_attributes_df
    pathway_entries = pathway_entries.replace("", np.nan).dropna(subset=["name"])

    # Map differentially expressed genes to the selected pathway
    data_df = merge_pathway_with_seq(pathway_entries, seq_df)

    # Select only connected nodes
    n_df, e_df = select_only_connected(data_df, pathway_edges)

    # Create a dictionary of the nodes
    n_df_dict = n_df.to_dict("index")

    # Add the name of the genes to the edge list
    e_df["source_name"] = e_df["source"].apply(lambda x: n_df_dict[x]["name_x"])
    e_df["target_name"] = e_df["target"].apply(lambda x: n_df_dict[x]["name_x"])

    g = nx.from_pandas_edgelist(
        e_df, "source", "target", edge_attr=True, create_using=nx.DiGraph()
    )

    nx.set_node_attributes(g, n_df.to_dict("index"))

    # Relabel the nodes indexing by gene names
    relabel_nodes(g)

    # Add attributes to the nodes such as degreee centrality
    # For directional graphs, this should be calculated as it's out-degree centrality
    deg_centrality = nx.out_degree_centrality(g)
    nx.set_node_attributes(g, deg_centrality, "degree_centrality")
    return g


def select_only_connected_from_name(merge_dataframe, edges_dataframe):

    _node_list = merge_dataframe["name"]

    edges_selected = edges_dataframe[
        (edges_dataframe["source"].isin(_node_list))
        & (edges_dataframe["target"].isin(_node_list))
    ]

    connect_nodes = list(edges_selected.source) + list(edges_selected.target)
    connect_nodes = list(set(connect_nodes))

    return (
        merge_dataframe[merge_dataframe["name"].isin(connect_nodes)].set_index("name"),
        edges_selected,
    )


# Formulate a final dataframe for visualization via Ensemble IDs
def combine_kegg_rna_seq_mouse(kegg_pth, rna_seq_pth):
    kegg_df = pd.read_csv(kegg_pth, index_col=0)
    rna_seq_df = pd.read_csv(rna_seq_pth, delimiter="\t")

    rna_seq_df["ens_gene"] = rna_seq_df["GeneID"].apply(lambda x: x.split(".")[0])
    df_mapped_to_kegg = pd.merge(
        kegg_df, rna_seq_df, how="inner", left_on="Ens", right_on="ens_gene"
    )
    df_mapped_to_kegg["fold_change"] = df_mapped_to_kegg["log2(FC)"].apply(
        lambda x: 2**x
    )

    return df_mapped_to_kegg


def generate_graph_list(root_path, kegg_path_list, species, sequence_path):
    """
    Generate a list of graphs from the provided KEGG path list and return pathway names.

    Returns:
        Tuple containing list of graphs and list of pathway names.
    """
    if species == "human":
        prefx = "hsa"
    elif species == "mouse":
        prefx = "mmu"
    elif species == "rat":
        prefx = "rno"
    else:
        raise ValueError("Only human, mouse and rat supported!")

    input_lst = [prefx + i for i in kegg_path_list]
    plst_name = [get_pathway_name(i) for i in input_lst]

    graph_list = []

    for name, pathway in zip(plst_name, input_lst):

        df_pth = os.path.join(root_path, f"{name}_{species}.csv")

        if os.path.exists(df_pth):
            if species == "mouse":
                mapped_df = combine_kegg_rna_seq_mouse(df_pth, sequence_path)
            else:
                mapped_df = combine_kegg_rna_seq(
                    df_pth, sequence_path, right_on="Gene name"
                )
        else:
            df = df_from_kegg_map(pathway)
            df["pathway"] = name
            rna_seq_df = pd.read_csv(sequence_path, delimiter="\t")

            if species == "mouse":
                rna_seq_df["ens_gene"] = rna_seq_df["GeneID"].apply(
                    lambda x: x.split(".")[0]
                )
                mapped_df = pd.merge(
                    df, rna_seq_df, how="inner", left_on="Ens", right_on="ens_gene"
                )
                mapped_df["fold_change"] = mapped_df["log2(FC)"].apply(lambda x: 2**x)
            else:
                mapped_df = pd.merge(
                    df, rna_seq_df, how="inner", left_on="name", right_on="Gene name"
                )
                mapped_df["fold_change"] = mapped_df["log2(FC)"].apply(lambda x: 2**x)

        g = get_graph_from_df(pathway, mapped_df)

        graph_list.append(g)

    return graph_list, plst_name


def save_node_importance(root_path, graphML):
    """
    Calculate and save node importance metrics (Degree Centrality, Eigenvector Centrality,
    Betweenness Centrality, Closeness Centrality, K-shell) to a specified path.

    Parameters
    ----------
    root_path : str
        Root path to save the node importance metrics.
    graphML : networkx.Graph
        GraphML object to calculate the node importance metrics from.

    Returns
    -------
    None
    """
    degDf = (
        pd.DataFrame.from_dict(
            nx.degree_centrality(graphML), orient="index", columns=["Degree Centrality"]
        )
        .reset_index()
        .rename(columns={"index": "geneID"})
    )
    degDf = degDf.sort_values(by="Degree Centrality", ascending=False)

    epiDf = (
        pd.DataFrame.from_dict(
            nx.eigenvector_centrality(graphML),
            orient="index",
            columns=["Eigenvector Centrality"],
        )
        .reset_index()
        .rename(columns={"index": "geneID"})
    )
    epiDf = epiDf.sort_values(by="Eigenvector Centrality", ascending=False)

    betDf = (
        pd.DataFrame.from_dict(
            nx.betweenness_centrality(graphML),
            orient="index",
            columns=["Betweenness Centrality"],
        )
        .reset_index()
        .rename(columns={"index": "geneID"})
    )
    betDf = betDf.sort_values(by="Betweenness Centrality", ascending=False)

    closDf = (
        pd.DataFrame.from_dict(
            nx.closeness_centrality(graphML),
            orient="index",
            columns=["Closeness Centrality"],
        )
        .reset_index()
        .rename(columns={"index": "geneID"})
    )
    closDf = closDf.sort_values(by="Closeness Centrality", ascending=False)

    shellDf = (
        pd.DataFrame.from_dict(
            nx.core_number(graphML), orient="index", columns=["K-shell"]
        )
        .reset_index()
        .rename(columns={"index": "geneID"})
    )
    shellDf = shellDf.sort_values(by="K-shell", ascending=False)

    node_importance_path = os.path.join(root_path, "important_node")

    if not os.path.exists(node_importance_path):
        os.makedirs(node_importance_path)

    degDf.to_csv(
        os.path.join(node_importance_path, "degree_centrality.csv"), index=False
    )
    epiDf.to_csv(
        os.path.join(node_importance_path, "eigenvector_centrality.csv"), index=False
    )
    betDf.to_csv(
        os.path.join(node_importance_path, "betweenness_centrality.csv"), index=False
    )
    closDf.to_csv(
        os.path.join(node_importance_path, "closeness_centrality.csv"), index=False
    )
    shellDf.to_csv(os.path.join(node_importance_path, "k_shell.csv"), index=False)
