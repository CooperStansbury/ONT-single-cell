import pandas as pd
import numpy as np
import os
import sys
import pyranges as pr

import anndata as ad
import scanpy as sc


def load_gene_table(fpath):
    """A function to load the gene table """
    gf = pr.read_gtf(fpath).df

    # subset just the major genes
    gf = gf[gf['Feature'] == 'gene']
    gf = gf[gf['gene_biotype'] == 'protein_coding'].reset_index(drop=True)
    print(f"{gf.shape=}")
    
    keep_cols = [
        'gene_id', 
        'gene_name',
        'Chromosome', 
        'Start',
        'End',
        'Strand',
    ]
    
    gf = gf[keep_cols]
    gf = gf[gf['gene_name'].notna()]
    gf = gf.drop_duplicates(subset='gene_name')
    return gf


def load_raw_counts(fpath):
    """A function to load the htseq-counts output"""
    df = pd.read_csv(fpath, sep='\t', low_memory=False)
    df = df.rename(columns={'Unnamed: 0' : 'gene_id'})
    return df


if __name__ == "__main__":
    counts_path = sys.argv[1]
    gene_path = sys.argv[2]
    out_path = sys.argv[3]

    # load the gene table
    gf = load_gene_table(gene_path)
    print(f"{gf.shape=}")

    # load the counts matrix
    df = load_raw_counts(counts_path)
    print(f"{df.shape=}")

    # filter out non-PT genes
    df = df[df['gene_id'].isin(gf['gene_id'])]
    df = df.set_index('gene_id')
    print(f"{df.shape=}")

    # make sure that the GTF has the same set
    gf = gf[gf['gene_id'].isin(df.index)]
    gf = gf.set_index('gene_id')
    print(f"{gf.shape=}")

    # TRANSPOSE for anndata
    df = df.T

    """Build andata object """
    adata = ad.AnnData(df.to_numpy())
    adata.obs_names = df.index
    adata.var_names = df.columns
    adata.var = gf

    # add gene counts as the first observation metric
    obs = df.sum(axis=1).reset_index(drop=False)
    obs.columns = ['cell_id', 'n_genes']
    obs = obs.set_index('cell_id')
    
    adata.obs = obs

    print(adata)
    print()
    
    # write the object to file
    adata.write(out_path)

    

    

    




    

   