import os
import sys
import pandas as pd
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
import scanpy as sc


def min_max(values):
    """Scales a series of values to the range [0, 1] using NumPy.

    Args:
        values: A Pandas Series or a NumPy array of numeric values.

    Returns:
        The scaled values as a Pandas Series or a NumPy array.
    """

    range_val = np.ptp(values)  # Peak-to-peak (max - min)
    if range_val == 0:
        return values
    
    scaled_values = (values - values.min()) / range_val
    return scaled_values


def anndata_stats(anndata, gene_list=None, label=None):
    """Calculates gene expression statistics for specified genes in an AnnData object.

    Args:
        anndata (anndata.AnnData): The AnnData object containing gene expression data.
        gene_list (list): A list of gene names to calculate statistics for.
        label (str, optional): An optional label to add to the results.

    Returns:
        pd.DataFrame: A DataFrame with gene statistics, including:
            - average_expression
            - median_expression
            - std_dev
            - num_nonzero_cells
            - percent_nonzero
            - mean_nonzero
            - median_nonzero
            - gene_name (if a label is provided)
    """
    
    if gene_list is None:
        gene_list = anndata.var_names

    # Validate gene list
    valid_gene_list = list(set(gene_list).intersection(anndata.var_names))
    if not valid_gene_list:
        raise ValueError("None of the provided genes are found in the AnnData object.")

    # Extract data and calculate statistics
    df = anndata[:, valid_gene_list].to_df()
    stats = pd.DataFrame({
        'average_expression': df.mean(),
        'median_expression': df.median(),
        'std_dev': df.std(ddof=0),
        'num_nonzero_cells': (df > 0).sum(),
        'percent_nonzero': 100 * (df != 0).mean(),
        'mean_nonzero': df[df > 0].mean(),
        'median_nonzero': df[df > 0].median(),
    })
    stats = stats.reset_index(names='gene_name')

    # Add optional label
    if label is not None:
        stats['label'] = label

    return stats


def calculate_gene_expression_stats(df):
    """
    Calculates aggregate gene expression statistics for a given DataFrame.

    Args:
        df (pd.DataFrame): A DataFrame with cells as rows and genes as columns.

    Returns:
        pd.DataFrame: A DataFrame with one row per gene, containing aggregated statistics.
    """

    # Aggregate statistics
    stats = pd.DataFrame({
        'average_expression': df.mean(),
        'median_expression': df.median(),
        'std_dev': df.std(ddof=0),  # Sample standard deviation
        'num_nonzero_cells': (df > 0).sum(),
        'percent_nonzero': 100 * (df != 0).mean()
    })

    # Round percent_nonzero to 2 decimal places
    stats['percent_nonzero'] = stats['percent_nonzero'].round(2)
    stats = stats.sort_values(by='average_expression', ascending=False)

    return stats


def top_n_de_genes(df, n=10, alpha=0.05, pct_nz_threshold=0.25, pct_nz_reference=0.90, values='logfoldchanges', gene_list=None):
    """
    Returns the top N differentially expressed genes for each group, ranked by absolute log fold change.

    Args:
        df (pd.DataFrame): DataFrame containing DE analysis results.
        n (int, optional): Number of top genes per group (default: 10).
        alpha (float, optional): Significance level (default: 0.05).
        pct_nz_threshold (float, optional): Min. proportion expressed in group (default: 0.25).
        pct_nz_reference (float, optional): Max. proportion expressed in reference (default: 0.90).
        values (str, optional): Pivot table values to fill (default: 'logfoldchanges')

    Returns:
        pd.DataFrame: A pivot table of log fold changes where rows are groups, columns are names, and values are logfoldchanges
    """

    # Filter by significance and expression percentage
    df_filtered = df.query(
        f"`pvals_adj` <= {alpha} and `pct_nz_group` >= {pct_nz_threshold} and `pct_nz_reference` <= {pct_nz_reference}"
    )
    
    if not gene_list is None:
        df_filtered = df_filtered[df_filtered['names'].isin(gene_list)]

    # Sort by log fold change within each group (descending order for top genes)
    df_sorted = df_filtered.sort_values(
        ['group', 'logfoldchanges'], ascending=[True, False]
    )

    # Group by and select top N rows
    top_genes = df_sorted.groupby('group').head(n)['names'].values

    table = df[df['names'].isin(top_genes)]

    table = pd.pivot_table(
        table,
        index='group',
        columns='names',
        values=values
    )
    
    table = table[top_genes]

    return table