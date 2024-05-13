import pandas as pd
import numpy as np
import os
import sys
import anndata as ad
import scanpy as sc

def basic_qc(adata, 
             min_genes=100, 
             min_cells=3, 
             target_sum=1e4,
             n_top_genes=2000):
    """
    Preprocesses single-cell RNA-seq data.

    Args:
        adata (sc.AnnData): The raw AnnData object
        min_genes (int, optional): Minimum number of genes expressed per cell. Defaults to 100.
        min_cells (int, optional): Minimum number of cells a gene is expressed in. Defaults to 3.
        target_sum (int, optional): Target sum for normalization. Defaults to 1e4.
        n_top_genes (int, optional): Number of highly variable genes to select. Defaults to 2000.

    Returns:
        sc.AnnData: The preprocessed AnnData object.
    """

    # Store raw counts
    adata.layers["raw_counts"] = adata.X.copy()

    # Filtering
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # Doublet detection
    sc.pp.scrublet(adata)
    
    # Saving count data
    adata.layers["filtered_counts"] = adata.X.copy()

    # Mitochondrial gene annotation
    adata.var['mt'] = adata.var['gene_name'].str.startswith('MT-')

    # QC metrics calculation
    sc.pp.calculate_qc_metrics(adata, 
                               qc_vars=['mt'], 
                               percent_top=None, 
                               log1p=False,
                               inplace=True)

    # Normalization and transformation
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # Highly variable gene selection 
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)

    return adata 


def default_embedding(
    adata, 
    n_pcs=50, 
    neighbors_k=15, 
    umap_min_dist=0.5, 
    n_neighbors=15
):
    """
    Performs PCA, calculates nearest neighbors, and computes UMAP embedding on an AnnData object.

    Args:
        adata: AnnData object containing the data to preprocess and embed.
        n_pcs: Number of principal components to compute. Default is 200.
        neighbors_k: Number of nearest neighbors to use for UMAP. Default is 15.
        umap_min_dist: Minimum distance between points in the UMAP embedding. Default is 0.5.
        n_neighbors: Number of neighbors to consider for nearest neighbor calculation. Default is 15.

    Returns:
        adata: The modified AnnData object with PCA, nearest neighbors, and UMAP results.
    """

    sc.tl.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    sc.tl.umap(adata, min_dist=umap_min_dist)

    return adata


def add_annotations(adata, dpath):
    """
    Maps gene names to Ensembl IDs and adds CSV annotations to an AnnData object.

    Args:
        adata: The AnnData object to which annotations will be added.
        dpath: The directory path containing CSV annotation files.

    Returns:
        The modified AnnData object with annotations added.
    """
    # add annother id column to var
    adata.var['ensembl_id'] = adata.var.index

    # Create gene name to Ensembl ID mapping
    gene_map = dict(zip(
        adata.var['gene_name'].astype(str).str.upper().values, 
        adata.var['ensembl_id'].values
    ))

    for filename in os.listdir(dpath):
        if not filename.endswith(".csv"):
            continue

        filepath = os.path.join(dpath, filename)
        df = pd.read_csv(filepath)

        # Ensure gene names are uppercase strings
        df['gene_name'] = df['gene_name'].astype(str).str.upper()

        # Filter to genes present in the gene_map
        df = df[df['gene_name'].isin(gene_map)]

        # Map gene names to Ensembl IDs
        df['ensembl_id'] = df['gene_name'].map(gene_map)

        # Add the dataframe to the AnnData object's uns attribute
        key_name = filename.replace(".csv", "")
        adata.uns[key_name] = df

    return adata



if __name__ == "__main__":
    anndata_path = sys.argv[1]
    out_path = sys.argv[2]
    annotation_directory = sys.argv[3]

    # load the data 
    adata = sc.read_h5ad(anndata_path)

    # process the data 
    adata = basic_qc(adata)
    
    # establish a deafault embedding
    adata = default_embedding(adata)

    # add gene annotations
    adata = add_annotations(adata, annotation_directory)
    
    # write the object to file
    adata.write(out_path)

    

    

    




    

   