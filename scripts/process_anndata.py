import pandas as pd
import numpy as np
import os
import sys
from types import SimpleNamespace  

import yaml
from pathlib import Path

import anndata as ad
import scanpy as sc



if __name__ == "__main__":
    anndata_path = sys.argv[1]
    config_path = sys.argv[2]
    out_path = sys.argv[3]

    print(f"{anndata_path=}")
    print(f"{config_path=}")
    print(f"{out_path=}")
    print()

    # load the config for processing params
    config = yaml.safe_load(Path(config_path).read_text())
    params = config['scanpy_params']

    print(params)
    print()

    # load the data 
    adata = sc.read_h5ad(anndata_path)
    adata.layers["counts"] = adata.X.copy()

    print(adata)
    print()
    
    # set up params
    min_genes = int(params['min_genes'])
    min_cells = int(params['min_cells'])
    target_sum = int(params['target_sum'])

    # perform simple filtering
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)



    # annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var['gene_name'].str.startswith('MT-')  

    # calc QC metrics
    sc.pp.calculate_qc_metrics(adata, 
                               qc_vars=['mt'], 
                               percent_top=None,
                               log1p=False, 
                               inplace=True)

    # normalize and transform
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    print(adata)
    
    # write the object to file
    adata.write(out_path)

    

    

    




    

   