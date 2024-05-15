import pandas as pd
import numpy as np
import os
import sys
import anndata as ad
import scanpy as sc

def float_to_string(float_num):
    """Converts a float between 0 and 1 to a zero-filled two-digit string."""
    return f"{int(float_num * 100):02d}"


if __name__ == "__main__":
    anndata_path = sys.argv[1]
    out_path = sys.argv[2]

    # load the data 
    adata = sc.read_h5ad(anndata_path)
    
    columns_to_drop = []
    for r in np.linspace(0, 1, num=100):
        r_val = float_to_string(r)
        key_added = f"r{r_val}"
        sc.tl.leiden(adata, 
                     resolution=r, 
                     key_added=key_added)
        
        columns_to_drop.append(key_added)
    
    # extract all clustering results and store in uns
    clusters = adata.obs[columns_to_drop].copy()
    adata.uns['clusters'] = clusters
        
    # drop the columns from adata.obs
    adata.obs = adata.obs.drop(columns_to_drop, axis=1)
    
    # remove metadata from uns
    for k in columns_to_drop:
        adata.uns.pop(k, None)
    
    # write the object to file
    adata.write(out_path)

    

    

    




    

   