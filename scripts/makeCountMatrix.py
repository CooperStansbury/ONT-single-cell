import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import scipy
    
    
def loadCounts(paths):
    """A function to load rna expression counts """
    res = []
    
    for fpath in paths:
        stem = Path(fpath).stem
        df = pd.read_csv(fpath, 
                         sep='\t', 
                         low_memory=False,
                         header=None)
        
        df.columns = ['geneId', stem]
        df = df.set_index('geneId')
        res.append(df)
        
    df = pd.concat(res, axis=1,)
    return df


def getGeneMap(gdf, protein_coding=True):
    """A function to get a hash table of gene ids
    and gene names. 
    
    NOTE: default is protien coding only
    """
    if protein_coding:
        gdf = gdf[gdf['gene_biotype'] == 'protein_coding']
        gdf = gdf.reset_index(drop=True)
    
    # drop unnamed genes
    gdf = gdf[gdf['gene_name'].notna()]
    gdf = gdf.reset_index(drop=True)

    # get the mapper
    geneMap = pd.Series(gdf.gene_name.values,index=gdf.gene_id).to_dict()
    return geneMap


def structureOutput(cdf, geneMap):
    """A function to structure the output, with a few
    assumptions """
    
    # drop the uncounted meta rows
    ndf = cdf[~cdf.index.str.contains("__")].copy()
    
    # map gene names
    ndf['geneName'] = ndf.index.map(geneMap)
    ndf = ndf[ndf['geneName'].notna()]
    
    # handle duplicated gene names with different IDs, 
    # for example haplotypes that aren't fully resolved
    ndf = cdf.groupby(ndf['geneName'], axis=0).sum() # this replaces the index
    return ndf

    



if __name__ == "__main__":
    outpath = sys.argv[1]
    genePath = sys.argv[2]
    inpaths = sys.argv[3:]
    
    # get gene names
    gdf = pd.read_csv(genePath, low_memory=False)
    
    # get gene names
    geneMap = getGeneMap(gdf, protein_coding=True)
    
    # load count information
    cdf = loadCounts(inpaths)
    
    # get the updated count matrix
    cdf = structureOutput(cdf, geneMap)
    cdf = cdf.reset_index(drop=False)
    cdf.to_csv(outpath, index=False)
    
