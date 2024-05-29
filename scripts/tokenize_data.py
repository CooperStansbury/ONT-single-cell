import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy.sparse as sp
import scanpy as sc
import anndata as an
from datasets import Dataset


KEEP_COLUMNS = [
    'cell_type',
    'n_counts',
    'dataset',
]


def load_gene_median_dict(gene_median_file):
    """
    Loads a gene median dictionary from a pickle file.

    Args:
        gene_median_file (str): Path to the pickle file containing the gene median dictionary.

    Returns:
        dict: A dictionary mapping gene IDs to their median expression values.
    """

    with open(gene_median_file, "rb") as f:
        gene_median_dict = pickle.load(f)

    return gene_median_dict


def load_gene_tokenization(token_dictionary_file):
    """
    Loads gene tokenization data from a pickle file.

    Args:
        token_dictionary_file (str): Path to the pickle file containing the gene-token dictionary.

    Returns:
        dict: Gene-token dictionary (Ensembl ID: token).
        list: List of all gene keys (Ensembl IDs).
        dict: Dictionary mapping gene keys to True (used for selecting genes later).
    """

    with open(token_dictionary_file, "rb") as f:
        gene_token_dict = pickle.load(f)

    gene_keys = list(gene_token_dict.keys())

    # Optimization: Pre-allocate the list for slight performance improvement
    genelist_dict = dict.fromkeys(gene_keys, True)

    return gene_token_dict, gene_keys, genelist_dict


def rank_genes(gene_vector, gene_tokens):
    """Ranks genes based on expression values in descending order.

    Args:
        gene_vector (numpy.ndarray): Array of gene expression values.
        gene_tokens (numpy.ndarray): Array of corresponding gene tokens.

    Returns:
        numpy.ndarray: Array of gene tokens sorted by descending expression value.
    """
    return gene_tokens[np.argsort(-gene_vector)]


def tokenize_anndata(adata, genelist_dict, gene_median_dict, 
                     chunk_size=10000, target_sum=10000, 
                     counts_column='n_counts', gene_id="ensemble_id"):
    """
    Tokenizes and ranks genes within an AnnData object, optimizing for memory efficiency.

    This function processes gene expression data in chunks, applies normalization, and ranks genes
    for each cell based on their expression levels. The resulting tokenized and ranked gene
    representations, along with cell metadata, are returned.

    Args:
        adata (AnnData): The AnnData object containing gene expression data.
        genelist_dict (dict): Dictionary mapping gene IDs to boolean values indicating relevance.
        gene_median_dict (dict): Dictionary mapping gene IDs to their median expression values.
        chunk_size (int, optional): Number of cells to process in each chunk (default: 1000).
        target_sum (int, optional): Target sum for count normalization (default: 10000).
        counts_column (str, optional): The column in `adata.obs` containing cell counts (default: 'n_counts').
        gene_id (str, optional): The column in `adata.var` containing gene IDs (default: 'ensembl_id').

    Returns:
        tuple: 
            - list: List of tokenized and ranked gene lists for each cell.
            - dict: Dictionary containing cell metadata (keys are metadata column names).
    """
    # Filter relevant miRNAs
    coding_miRNA_mask = np.array([genelist_dict.get(i, False) for i in adata.var[gene_id]])
    coding_miRNA_loc = np.where(coding_miRNA_mask)[0]

    # Extract miRNA information
    coding_miRNA_ids = adata.var[gene_id][coding_miRNA_loc]
    norm_factor_vector = np.array([gene_median_dict[i] for i in coding_miRNA_ids])
    coding_miRNA_tokens = np.array([gene_token_dict[i] for i in coding_miRNA_ids])

    tokenized_cells = []
    file_cell_metadata = {k: [] for k in adata.obs.columns}  # Initialize metadata dict

    # Process in chunks for memory efficiency
    for chunk_start in range(0, adata.shape[0], chunk_size):
        chunk_end = chunk_start + chunk_size
        adata_chunk = adata[chunk_start:chunk_end, coding_miRNA_loc]

        # Normalize counts
        n_counts = adata_chunk.obs[counts_column].values[:, None]
        X_norm = adata_chunk.X / n_counts * target_sum / norm_factor_vector
        X_norm = sp.csr_matrix(X_norm)  

        # Tokenize and rank genes for each cell in chunk
        for i in range(X_norm.shape[0]):
            ranks = rank_genes(X_norm[i].data, coding_miRNA_tokens[X_norm[i].indices])
            ranks = list(ranks[~np.isnan(ranks)].astype(int))

            tokenized_cells.append(ranks)

        # Update metadata
        for k in adata.obs.columns:
            file_cell_metadata[k].extend(adata_chunk.obs[k].tolist())

    return tokenized_cells, file_cell_metadata


def create_dataset(tokenized_cells, cell_metadata, gene_token_dict, model_input_size=2048):
    """
    Creates a Hugging Face Dataset from tokenized cells and associated metadata.

    Args:
        tokenized_cells (list): List of tokenized cell representations (lists of tokens).
        cell_metadata (dict, optional): Dictionary containing additional cell metadata.
        model_input_size (int): The maximum input size for the model.
        gene_token_dict (dict): Dictionary mapping genes to their tokens.

    Returns:
        datasets.Dataset: The processed Hugging Face dataset.
    """
    
    # Merge cell metadata into the dataset dictionary if provided
    dataset_dict = {
        "input_ids": tokenized_cells,
        **cell_metadata 
    }
    
    output_dataset = Dataset.from_dict(dataset_dict)

    def format_cell_features(example):
        example["input_ids"] = example["input_ids"][0 : model_input_size] # truncate
        example["length"] = len(example["input_ids"])  # Add length for convenience
        return example

    nproc = multiprocessing.cpu_count()
    return output_dataset.map(format_cell_features, num_proc=nproc)  # Return mapped dataset


def save_hf_dataset(dataset: Dataset, output_path: str, overwrite=True):
    """
    Saves a Hugging Face Dataset to disk at a specified file path.

    This function serializes a Hugging Face `Dataset` object and saves it to disk in the Arrow format.

    Args:
        dataset (Dataset): The Hugging Face `Dataset` object to be saved.
        output_path (str): The full file path (including the filename) where the dataset will be saved. 
        overwrite (bool, optional): If `True`, an existing dataset at `output_path` will be overwritten. 
                                   If `False` and the file exists, a `FileExistsError` is raised (default: True).

    Raises:
        TypeError: If `dataset` is not a Hugging Face `Dataset` instance.
        FileExistsError: If `output_path` points to an existing file and `overwrite` is False.
    """

    if not isinstance(dataset, Dataset):
        raise TypeError("The provided dataset is not a Hugging Face Dataset.")

    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError(
            f"Dataset '{output_path}' already exists. Set `overwrite=True` to overwrite."
        )
    dataset.save_to_disk(output_path)

        
if __name__ == "__main__":
    adata_path = sys.argv[1]
    output_path = sys.argv[2]
    gene_median = sys.argv[3]
    token_dict = sys.argv[4]
    
    # load the data 
    gene_token_dict, gene_keys, genelist_dict = load_gene_tokenization(token_dict)
    gene_median_dict = load_gene_median_dict(gene_median)
    
    # load and restructure
    adata = sc.read_h5ad(adata_path)
    
    # drop most metadata columns
    adata.obs = adata.obs[KEEP_COLUMNS]

    # tokenize
    tokenized_cells, cell_metadata = tokenize_anndata(adata, 
                                                      genelist_dict, 
                                                      gene_median_dict)
    
    # Merge cell metadata into the dataset dictionary
    dataset_dict = {
        "input_ids": tokenized_cells,
        **cell_metadata 
    }
    
    output_dataset = Dataset.from_dict(dataset_dict)

    def format_cell_features(example):
        example["input_ids"] = example["input_ids"][0 : model_input_size] # truncate
        example["length"] = len(example["input_ids"])  # Add length for convenience
        return example

    dataset = output_dataset.map(format_cell_features, num_proc=16) 
    
    # store output
    save_hf_dataset(dataset, output_path, overwrite=True)
    
    
    
    
    
    
    