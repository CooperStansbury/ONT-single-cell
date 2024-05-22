import pandas as pd
import numpy as np
import os
import sys
import anndata as ad
import scanpy as sc
from scipy import sparse

def basic_qc(adata, 
             min_genes=100, 
             min_cells=3, 
             target_sum=1e4):
    """
    Preprocesses single-cell RNA-seq data.

    Args:
        adata (sc.AnnData): The raw AnnData object
        min_genes (int, optional): Minimum number of genes expressed per cell. Defaults to 100.
        min_cells (int, optional): Minimum number of cells a gene is expressed in. Defaults to 3.
        target_sum (int, optional): Target sum for normalization. Defaults to 1e4.

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
    return adata 


def get_embeddings(adata):
    """
    Performs standard preprocessing and embedding steps for single-cell data.

    This function applies PCA, computes nearest neighbors, and generates both UMAP 
    and t-SNE embeddings on an AnnData object.

    Args:
        adata (AnnData): The AnnData object containing the data.

    Returns:
        AnnData: The modified AnnData object with the calculated embeddings.
    """

    # PCA (Principal Component Analysis)
    sc.tl.pca(adata, svd_solver='arpack')  # 'arpack' is often more efficient for large datasets

    # Nearest Neighbors
    sc.pp.neighbors(adata)

    # UMAP (Uniform Manifold Approximation and Projection)
    sc.tl.umap(adata)

    # t-SNE (t-Distributed Stochastic Neighbor Embedding) 
    sc.tl.tsne(adata)  # t-SNE is now always computed

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
    
    print("------------------------ Raw Matrix ------------------------")
    sc.logging.print_memory_usage()

    # process the data 
    adata = basic_qc(adata)
    
    # Store raw counts
    adata.layers["log_norm"] = adata.X.copy()
    
    print("------------------------ Post QC ------------------------")
    sc.logging.print_memory_usage()
    
    """BATCH CORRECTION"""
    sc.pp.combat(adata, key='dataset')

    # handle negatives
    adata.X = sparse.csr_matrix(np.where(adata.X < 0, 0, adata.X))
    
    print("------------------------ Post Combat ------------------------")
    sc.logging.print_memory_usage()
    
    # Highly variable gene selection 
    sc.pp.highly_variable_genes(adata)
    
    # establish a deafault embeddings
    adata = get_embeddings(adata)

    # add gene annotations
    adata = add_annotations(adata, annotation_directory)
    
    print("------------------------ Post Annotation ------------------------")
    sc.logging.print_memory_usage()
    
    # write the object to file
    adata.write(out_path)

    

    

    




    

   