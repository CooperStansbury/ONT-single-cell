import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse
import anndata as an
import gc 

def clean_adata(adata, additional_obs=None, additional_var=None, additional_uns=None):
    """
    Cleans an AnnData object by removing specified data.

    Args:
        adata: The AnnData object to clean.
        additional_obs: List of additional obs keys to remove (optional).
        additional_var: List of additional var keys to remove (optional).
        additional_uns: List of additional uns keys to remove (optional).

    Returns:
        The cleaned AnnData object.
    """

    # Data to remove
    data_to_remove = {
        'obsp': ['connectivities', 'distances'],
        'obsm': ['X_pca', 'X_umap', 'X_tsne'],
        'uns': [
            'tsne', 'pca', 'umap', 'scrublet', 'neighbors', 'HSC_vs_FB', 'HSC_vs_FB_pure', 
            'fb_vs_hsc_up', 'hsc_v_fib_up', 'hvg', 'log1p', 'scenic_transcription_factors', 
            'tabula_sapiens_deg'
        ],
        'varm': ['PCs'],
        'var': [
            'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts',
            'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 
            'ensembl_id'
        ],
        'obs': [
            'n_genes', 'doublet_score', 'predicted_doublet', 'n_genes_by_counts', 'total_counts',
            'total_counts_mt', 'pct_counts_mt', 'Barcode', 'Library', 'method', 'donor', 
            'anatomical_information', 'n_counts_UMIs', 'cell_ontology_class', 'free_annotation',
            'manually_annotated', 'compartment', 'gender', 'celltype', 'record_id', 'cell_id',
            'nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC', 'nCount_SCT', 
            'nFeature_SCT', 'SCT.weight', 'ATAC.weight', 'seurat_clusters', 'STD.CellType', 
            'STD_Cat', 'STD_Cat2', 'Sample', 'HLF', 'CRHBP', 'CD34', 'MitoCoverage', 'ClonalGroup', 
            'Sig.HSC1', 'Sig.Prog1', 'Sig.EarlyE1', 'Sig.LateE1', 'Sig.ProMono1', 'Sig.Mono1', 
            'Sig.ncMono1', 'Sig.cDC1', 'Sig.pDC1', 'Sig.ProB1', 'Sig.PreB1', 'Sig.B1', 'Sig.Plasma1',
            'Sig.T1', 'Sig.CTL1', 'Sig.NK1', 'meanCov', 'ClonalGroup.Prob', 'wsnn_res.0.8', 
            'Origin.Seurat'
        ],
        'layers': [
            'log_norm',
        ]
        
    }

    # Add additional data to remove if provided
    if additional_obs:
        data_to_remove['obs'].extend(additional_obs)
    if additional_var:
        data_to_remove['var'].extend(additional_var)
    if additional_uns:
        data_to_remove['uns'].extend(additional_uns)

    # Remove the data
    for key, items in data_to_remove.items():
        for item in items:
            try:
                del getattr(adata, key)[item]
            except KeyError:
                pass  # Skip if the key doesn't exist
    gc.collect()
    return adata


        
if __name__ == "__main__":
    adata_path = sys.argv[1]
    output_path = sys.argv[2]
  
    # load the "heavy" ann data object
    adata = sc.read_h5ad(adata_path)
    
    print("------------------------ Raw Matrix ------------------------")
    sc.logging.print_memory_usage()
    
    adata = clean_adata(adata)
    print("------------------------ Trimmed Matrix ------------------------")
    sc.logging.print_memory_usage()
    
    # drop some cell types
    keep = [
        'iHSC',
        'LinNegCD34lowCD164high',
        'HSC',
        'LinNegCD34PosCD164Pos',
        'MPP',
        'MLP',
        'FB',
        'MKP',
        'Refined.HSC',
        'LMPP',
    ]

    adata = adata[adata.obs['cell_type'].isin(keep), :]
    gc.collect()
    
    # drop some datasets
    data = [
        'pellin', # pelling 2019
        'iHSC', # our iHSC data
        'old1_BMMC_HSPC', # weng data
        'old2_BMMC_HSPC', # weng data
        'young2_HSC', # weng data
        'tabula_sapiens', # fibroblast
    ]

    adata = adata[adata.obs['dataset'].isin(data), :]
    gc.collect()
    
    # drop HSC from tabula
    mask = (adata.obs['cell_type'] == 'HSC') & (adata.obs['dataset'] == 'tabula_sapiens')
    adata = adata[~mask, :]
    gc.collect()
    
    print("------------------------ Subsetted Matrix ------------------------")
    sc.logging.print_memory_usage()
    
    # re-process and combat
    adata.X = adata.layers['raw_counts']
    
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=100)
    
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers["log_norm"] = adata.X.copy()
    
    # batch correction
    sc.pp.combat(adata, key='dataset')

    # handle negatives
    adata.X = sparse.csr_matrix(np.where(adata.X < 0, 0, adata.X))
    
    # make matrix sparse, reduce the precision
    adata.X = sparse.csr_matrix(adata.X.astype('float32'))
    
    print("------------------------ Post Downcasting ------------------------")
    sc.logging.print_memory_usage()
    
    # output the results
    adata.copy().write(output_path)
    
    
    
    
    
    
    
    

 
    
    
    
    