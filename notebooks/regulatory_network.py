import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


"""SCENIC FUNCTIONS"""

def load_scenic_database(filepath="/nfs/turbo/umms-indikar/shared/projects/HSC/data/scenic_resources/500bp_up_100bp_down_B.csv"):
    """Loads the Scenic database from a CSV file.

    Args:
        filepath (str, optional): The path to the CSV file. 

    Returns:
        pandas.DataFrame: The Scenic database with 'gene_name' as the index.
    """

    sdf = pd.read_csv(filepath)
    sdf = sdf.rename(columns={"Unnamed: 0": "gene_name"})
    sdf = sdf.set_index("gene_name")
    return sdf


def get_regulators(sdf, query, num_results=10):
    """Retrieves the top regulators of a gene or set of genes in a Scenic dataframe.

    Args:
        sdf (pd.DataFrame): The Scenic dataframe.
        query (str or list): The gene(s) of interest to find regulators for.
        num_results (int, optional): The number of top regulators to return. Defaults to 10.

    Returns:
        pd.Index: An index containing the names of the top regulators.
    """
    return sdf.loc[query].nlargest(num_results).index.to_list()


def get_targets(sdf, query, num_results):
    """Extract the top binding targets of a given transcription factor.

    Args:
        sdf (pd.DataFrame): Scenic database DataFrame.
        query (str): Column name representing the transcription factor.
        num_results (int): Number of top targets to return.

    Returns:
        list: Top binding targets.
    """
    # Sort values in descending order and get the top 'num_results' indices
    top_target_indices = sdf[query].nlargest(num_results).index
    return top_target_indices.tolist()


def get_targets_quantile(sdf, query, quantile):
    """Extract binding targets of a given transcription factor, filtered by quantile.

    Args:
        sdf (pd.DataFrame): Scenic database DataFrame.
        query (str): Column name representing the transcription factor.
        quantile (float): Quantile threshold (0.0 to 1.0) for filtering.

    Returns:
        list: Binding targets exceeding the quantile threshold.
    """
    query_values = sdf[query]
    threshold = query_values.quantile(quantile)
    filtered_targets = query_values[query_values > threshold]
    return filtered_targets.index.tolist() 


def create_scenic_network(gene_expression_df, scenic_tf_matrix_df):
    """
    Creates a directed gene regulatory network based on transcription factor (TF) activity 
    from a SCENIC analysis.

    Args:
        gene_expression_df (pd.DataFrame): DataFrame with gene expression data.
        scenic_tf_matrix_df (pd.DataFrame): DataFrame with SCENIC TF-target gene interactions.

    Returns:
        nx.DiGraph: Directed network where nodes are genes and edges represent TF-target gene interactions.
    """
    
    # Filter TFs from the gene expression data
    active_tfs = scenic_tf_matrix_df.columns.intersection(gene_expression_df["gene_name"].unique())
    gene_expression_df['is_tf'] = gene_expression_df['gene_name'].isin(active_tfs)

    # Subset SCENIC matrix to active TFs and expressed genes, reshape to edgelist format
    sdf = scenic_tf_matrix_df.loc[scenic_tf_matrix_df.index.isin(gene_expression_df["gene_name"]), active_tfs]
    sdf = sdf.stack().reset_index()
    sdf.columns = ["gene_name", "TF", "weight"]

    # Create directed graph using NetworkX
    G = nx.from_pandas_edgelist(
        sdf, source="TF", target="gene_name", edge_attr="weight", create_using=nx.DiGraph()
    )

    # Map gene expression data to the network nodes
    nx.set_node_attributes(G, gene_expression_df.set_index("gene_name").to_dict(orient="index"))
    
    return G


""" GENERIC PLOTTING """

def visualize_network(G, **kwargs):
    """
    Visualizes a NetworkX graph G using Kamada-Kawai layout. Accepts additional
    keyword arguments to pass directly to `nx.draw_networkx`.

    Args:
        G (nx.Graph): The NetworkX graph to visualize.
        **kwargs: Additional keyword arguments to pass to `nx.draw_networkx`.
    """

    # Kamada-Kawai layout (assumed default)
    pos = nx.kamada_kawai_layout(G)

    # Extract edges
    edges = list(G.edges())

    # Determine node colors based on 'is_tf' attribute (or provide your own logic)
    colors = ["lightcoral" if G.nodes[n].get("is_tf", False) else "lightskyblue" for n in G.nodes]

    # Edge weights based on "weight" attribute (or provide your own logic)
    edge_weights = np.array(list(nx.get_edge_attributes(G, "weight").values()))
    edge_weights = np.log1p(edge_weights ** 1.1) + 0.1

    # Default parameters (can be overridden by kwargs)
    default_params = {
        "node_color": colors,
        "alpha": 1,
        "edgecolors": "k",
        "node_size": 800,
        "linewidths": 1.0,
        "with_labels": True,
        "font_color": "k",
        "font_size": 5,
        "font_family": "sans-serif",
        "font_weight": "bold",
        "horizontalalignment": "center",
        "verticalalignment": "center",
        "edgelist": edges,
        "width": edge_weights,
        "edge_color": "k",
        "arrows": True,
        "arrowsize": 10
    }

    # Update parameters with user-provided kwargs
    params = default_params | kwargs  # Merge dictionaries (Python 3.9+ syntax)

    # Draw the graph
    nx.draw_networkx(G, pos, **params)