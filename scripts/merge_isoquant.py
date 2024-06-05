import pandas as pd
import numpy as np
import os
import sys

def combine_gene_counts(file_list):
    """Combines gene count data from multiple tab-separated files.

    Args:
        file_list (list): A list of file paths to the tab-separated files.

    Returns:
        pd.DataFrame: A DataFrame with gene IDs as the index, columns representing each file's counts,
                      and a 'total_counts' column summarizing counts across all files.
    """
    dfs = []

    for fpath in file_list:
        run_id = os.path.basename(fpath).split(".")[0]  # Extract run ID from filename
        tmp_df = pd.read_csv(
            fpath, 
            sep="\t", 
            names=['gene_id', run_id],  # Use run ID as column name
            comment="#",               # Ignore lines starting with "#"
            index_col='gene_id'        # Set 'gene_id' as the index
        )
        dfs.append(tmp_df)

    # Concatenate DataFrames and Calculate Total Counts
    df = pd.concat(dfs, axis=1)   
    df['total_counts'] = df.sum(axis=1)  # Sum across columns to get total counts

    return df

        
if __name__ == "__main__":
    outpath = sys.argv[1]
    file_list = sys.argv[2:]
    
    df = combine_gene_counts(file_list)
    df = df.reset_index()
    df.to_csv(outpath, index=False)
    
    


   