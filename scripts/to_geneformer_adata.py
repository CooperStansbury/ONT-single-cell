import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as an


def load_gene_map(gene_table_path, gene_type="protein_coding", key="gene_name"):
    """
    Loads a gene mapping table from a CSV file, filtering by gene type and creating a dictionary.

    This function reads gene information from a CSV file, filters for the specified gene type 
    (defaulting to protein-coding genes), removes duplicates and missing values, and then 
    constructs a dictionary mapping the specified column (`key`) to its corresponding values.
    
    The default behavior is to map gene names ('gene_name') to gene IDs ('gene_id'), but this can be
    customized by changing the `key` argument.

    Args:
        gene_table_path (str): The path to the CSV file containing the gene table.
        gene_type (str, optional): The type of gene to filter for (default: "protein_coding").
        key (str, optional): The column to use as keys in the dictionary (default: "gene_name").

    Returns:
        dict: A dictionary mapping gene keys (as specified by `key`) to the corresponding gene values.
    """
    usecols = ['gene_id', 'gene_name', 'gene_biotype']
    df = pd.read_csv(gene_table_path, usecols=usecols)

    # Filter and clean data in a single chain
    df = (
        df.drop_duplicates()
        .query("gene_biotype == @gene_type")
        .dropna(subset=['gene_name', 'gene_id'])
    )
    
    if key == 'gene_name':
        return dict(zip(df['gene_name'], df['gene_id']))
    elif key == 'gene_id':
        return dict(zip(df['gene_id'], df['gene_name']))
            

def skeletonize(adata, gene_map, gene_identifier, gene_column_type, gene_index, counts_layer):
    """
    Creates a simplified AnnData object by filtering genes and mapping identifiers.

    This function takes an AnnData object and a gene mapping dictionary. It filters the genes in 
    the AnnData based on the provided mapping, optionally converting between gene names and Ensembl IDs.
    The resulting AnnData object contains only the relevant genes and their associated count data,
    along with the original observation metadata.

    Args:
        adata (anndata.AnnData): The input AnnData object containing gene expression data.
        gene_map (dict): A dictionary mapping gene identifiers (keys) to desired values.
                            Keys can be either gene names or Ensembl IDs, depending on `gene_column_type`.
        gene_identifier (str): The column name in `adata.var` containing the gene identifiers 
                            (either 'gene_name' or 'ensemble_id').
        gene_column_type (str): Indicates the type of gene identifier in the `gene_identifier` column:
                            'gene_name' or 'ensemble_id'.
        gene_index (str): The column name to use as the gene index in the output AnnData's `.var` attribute.
                            Typically 'gene_name' or 'ensemble_id'.
        counts_layer (str): The layer name in `adata` to use for count data (default is 'counts').

    Returns:
        anndata.AnnData: A new AnnData object with:
            - `.X`: Count data filtered to the relevant genes.
            - `.obs`: Copied from the original `adata`.
            - `.var`: Contains the `gene_index` as the index, and the mapped gene identifiers.
            - `.obs['n_counts']`: Added (or preserved) to indicate the total counts per cell.

    Raises:
        ValueError: If an invalid `gene_column_type` is provided.
    """
    
    if gene_column_type == 'gene_name':
        var = pd.DataFrame(adata.var[gene_identifier].copy())
        var.columns = ['gene_name']
        var['ensemble_id'] = var['gene_name'].map(gene_map)
        var = var[var['ensemble_id'].notna()]
        gene_idx = var.index.astype(str) # use the existing index
    
    elif gene_column_type == 'ensemble_id':
        var = pd.DataFrame(adata.var[gene_identifier].copy())
        var.columns = ['ensemble_id']
        var['gene_name'] = var['ensemble_id'].map(gene_map)
        var = var[var['gene_name'].notna()]
        gene_idx = var.index.astype(str) # use the existing index
    else:
        raise ValueError('gene_column_type must be one of: `ensemble_id` or `gene_name`')
        
    X = adata[:, gene_idx].layers[counts_layer].copy()
    ndata = an.AnnData(X)
    ndata.obs = adata.obs.copy()
    ndata.obs_names = adata.obs_names
    
    if not 'n_counts' in ndata.obs.columns:
        ndata.obs['n_counts'] = X.sum(axis=1)
    
    ndata.var_names = var[gene_index].astype(str).values
    ndata.var = var.set_index(gene_index)
    ndata.var_names_make_unique()
    
    return ndata

        
if __name__ == "__main__":
    adata_path = sys.argv[1]
    gene_table_path = sys.argv[2]
    gene_table_key = sys.argv[3]
    gene_identifier = sys.argv[4] 
    gene_column_type = sys.argv[5] 
    gene_index = sys.argv[6] 
    counts_layer = sys.argv[7] 
    output_path = sys.argv[8]
    
    # load the data
    adata = sc.read_h5ad(adata_path)
    
    # load the gene map
    gene_map = load_gene_map(gene_table_path, key=gene_table_key)
    
    # create the new data structure
    ndata = skeletonize(adata, gene_map,
                       gene_identifier=gene_identifier, 
                       gene_column_type=gene_column_type, 
                       gene_index=gene_index,
                       counts_layer=counts_layer)
    
    # save the data
    ndata.write(output_path)
    
    

 
    
    
    
    