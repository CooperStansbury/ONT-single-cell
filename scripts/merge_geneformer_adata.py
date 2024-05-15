import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as an

CELL_TYPE_COLUMNS = ['celltype', 'STD.CellType', 'dataset']


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
            

def consolidate_cell_type(
    df: pd.DataFrame, celltype_columns: list, fill_value: str = 'iHSC'
) -> pd.DataFrame:
    """Consolidates cell type information from multiple columns into a single 'cell_type' column.

    Args:
        df: The DataFrame containing the cell type data.
        celltype_columns: List of column names to check for cell type values.
        fill_value (str, optional): The value to fill in 'cell_type' if none of the columns in `celltype_columns` exist in the DataFrame (default: 'iHSC').

    Returns:
        The DataFrame with a new 'cell_type' column containing consolidated values or `fill_value` if no overlap is found.
    """
    
    # Filter columns that exist in the DataFrame
    valid_columns = [col for col in celltype_columns if col in df.columns]

    if valid_columns:  
        # Find the first non-null value from the specified columns, row-wise
        df['cell_type'] = df[valid_columns].bfill(axis=1).iloc[:, 0]
    else:
        # If no valid columns, fill 'cell_type' with the provided fill_value
        df['cell_type'] = fill_value
    
    # Replace NaN with None to ensure the new column has the proper data type
    df['cell_type'] = df['cell_type'].astype('object')
    df['cell_type'] = df['cell_type'].where(pd.notnull(df['cell_type']), None)
    
    return df

        
if __name__ == "__main__":
    gene_table_path = sys.argv[1]
    output_path = sys.argv[2]
    file_list = sys.argv[3:]
    
    # load the gene map
    gene_map = load_gene_map(gene_table_path)
    
    # collect the data
    adatas = {}

    for fpath in file_list:
        base_name = os.path.basename(fpath)
        file_id = base_name.replace(".anndata.h5ad", "")
        tmp = sc.read_h5ad(fpath)
        tmp.var_names_make_unique()
        
        # handle issues with categoricals
        not_na = (tmp.var['ensemble_id'].notna()) & (tmp.var['ensemble_id'].isin(list(gene_map.values())))
        tmp = tmp[:, not_na].copy()
        tmp.var['ensemble_id'] = tmp.var['ensemble_id'].astype(str)
        
        tmp.obs = consolidate_cell_type(tmp.obs, CELL_TYPE_COLUMNS)
        adatas[file_id] = tmp
        
    # merge the data
    adata = an.concat(adatas, 
                      label='dataset', 
                      index_unique="_", 
                      join="outer")

    # set up the variable names
    adata.var['gene_name'] = adata.var.index
    adata.var_names_make_unique()
    adata.var['ensemble_id'] = adata.var['gene_name'].map(gene_map)
    
    # drop non-ensemble genes
    not_na = adata.var['ensemble_id'].notna()
    adata = adata[:, not_na]
    
    # Convert all categorical columns to strings
    for col in adata.obs.select_dtypes(include='category'):
        adata.obs[col] = adata.obs[col].astype(str)
    
    # output the results
    adata.write(output_path)
    
    
    
    
    
    
    
    

 
    
    
    
    