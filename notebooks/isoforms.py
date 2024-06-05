import os
import numpy as np
import pandas as pd
import gget
import glob
import networkx as nx
from scipy.stats import entropy
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


def load_isoquant(
    fpath: str,
    biotype: str = "protein_coding",
    counts_col: str = "gene_CPM",
    min_transcript_counts: int = 1,
    min_isoforms_per_gene: int = 2,
    exclude_prefixes: list[str] = ["MT", "RP"],
) -> pd.DataFrame:
    """
    Loads IsoQuant data from a CSV file, filters, and transforms it for further analysis.

    Args:
        fpath (str): Path to the IsoQuant CSV file.
        biotype (str): Keep only transcripts of this biotype (default: "protein_coding").
        counts_col (str): Name of the column containing gene expression values (default: "gene_CPM").
        min_transcript_counts (int): Minimum transcript count for inclusion (default: 1).
        min_isoforms_per_gene (int): Minimum number of isoforms per gene to keep (default: 2).
        exclude_prefixes (list[str]): Exclude genes whose names start with these prefixes.

    Returns:
        pd.DataFrame: A DataFrame containing the filtered and processed IsoQuant data, with a new 'tid' 
                      column representing transcript IDs within each gene.
    """
    # load 
    df = pd.read_csv(fpath, low_memory=False)
    
    # Filtering
    df = df.dropna(subset=['gene_name', 'transcript_name'])
    df = df[df['transcript_biotype'] == biotype]
    for prefix in exclude_prefixes:
        df = df[~df['gene_name'].str.startswith(prefix)]

    # Transformations
    df['log_CPM'] = np.log1p(df["gene_CPM"])

    # Filtering on Expression and Isoform Count
    df = df[df['transcript_counts'] > min_transcript_counts]
    df['n_isoforms'] = df.groupby('gene_name')['transcript_name'].transform('nunique')
    df = df[df['n_isoforms'] >= min_isoforms_per_gene]

    # Sorting and Transcript ID Assignment
    df = df.sort_values(by=['gene_name', 'transcript_name'])
    df['tid'] = df.groupby('gene_name').cumcount() + 1
    
    # Create the pattern to remove
    df['pattern'] = df['gene_name'].astype(str) + '-'
    # Remove the pattern from the transcript name
    df['short_name'] = df.apply(lambda row: row['transcript_name'].replace(row['pattern'], '', 1), axis=1)
    df = df.drop(columns='pattern')

    return df


def extract_gene_data(pdf):
    """Extracts gene-level information from a pandas DataFrame.

    Args:
        pdf: A pandas DataFrame containing gene-level and transcript-level data.

    Returns:
        A pandas DataFrame with only unique gene-level information.
    """

    cols = ['gene_name', 'gene_biotype', 'gene_counts', 
            'gene_CPM', 'log_CPM', 'n_isoforms', 'entropy']
    genes_only = pdf[cols].drop_duplicates()
    return genes_only


def calculate_tx_percent(pdf, genes=None, count_threshold=1000):
    """
    Calculates the percentage of each transcript's expression relative to its gene's total expression,
    and returns a DataFrame with the calculated percentages.

    Args:
        pdf (pd.DataFrame): DataFrame with gene and transcript information,
                             including columns 'gene_name', 'transcript_counts', and 'gene_counts'.
        genes (list, optional): List of genes to include. If None, uses all unique genes.
        count_threshold (int, optional): Minimum gene count threshold. Defaults to 1000.

    Returns:
        pd.DataFrame: A DataFrame containing the original columns and a new 'tx_percent' column.
    """

    if genes is None:
        genes = pdf['gene_name'].unique()

    # Filter and calculate tx_percent in one step
    filtered_pdf = (
        pdf.loc[
            (pdf['gene_name'].isin(genes)) & (pdf['gene_counts'] > count_threshold)
        ]
        .assign(
            tx_percent=lambda df: df['transcript_counts']
            / df.groupby('gene_name')['transcript_counts'].transform('sum')
        )
    )

    # Check if filtering resulted in any valid data
    if filtered_pdf.empty:
        raise ValueError("No transcripts found matching the filtering criteria.")

    return filtered_pdf


def pivot_isoforms(pdf, translate_column='short_name', genes=None, count_threshold=1000, drop_zero=True):
    """
    Creates a pivot table of transcript percentages per gene, based on filtering criteria.

    This function first filters the input DataFrame based on specified genes and a count threshold.
    Then, it calculates transcript percentages and creates a pivot table with genes as rows, 
    transcript labels as columns, and percentages as values.

    Args:
        pdf (pd.DataFrame): DataFrame with gene and transcript information.
        translate_column (str): Column to use for transcript labels (default: 'short_name').
        genes (list): List of genes to include (default: None, includes all unique genes).
        count_threshold (int): Minimum gene count threshold (default: 1000).
        drop_zero (bool): Whether to drop columns with all zero values (default: True).

    Returns:
        tuple: Tuple containing:
            * table (pd.DataFrame): Pivot table with transcript percentages per gene.
            * filtered_pdf (pd.DataFrame): Filtered DataFrame used to create the pivot table.
    """
    # Apply Filters Directly to DataFrame
    filtered_pdf = calculate_tx_percent(pdf, genes=genes, 
                                        count_threshold=count_threshold)
    
    # Clearer Pivot Table Creation using 'fillna'
    table = filtered_pdf.pivot_table(
        index='gene_name',
        columns=translate_column,
        values='tx_percent'
    ).fillna(0.0)  # 0.0 is more explicit for percentage values
    

    if drop_zero:
        zero_sum_cols = table.columns[(table.sum(axis=0) == 0)]
        table = table.drop(columns=zero_sum_cols)


    return table, filtered_pdf


def calculate_gene_entropy(df):
    """
    Calculates the entropy of transcript expression for each gene and adds it as a new column to the DataFrame.

    Args:
        df (pd.DataFrame): DataFrame with columns 'gene_name' and 'tx_percent'.

    Returns:
        pd.DataFrame: The original DataFrame with an added 'entropy' column.
    """

    # Validate required columns
    required_cols = ['gene_name', 'tx_percent']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"DataFrame must contain the following columns: {required_cols}")

    def gene_entropy(group):
        """Calculates the entropy for a group of transcripts belonging to a single gene."""
        return entropy(group[group > 0], base=2)

    # Calculate entropy for each gene group and merge back to the original DataFrame
    entropy_series = df.groupby("gene_name")["tx_percent"].agg(gene_entropy)
    df = df.merge(entropy_series.rename("entropy"), on="gene_name", how="left")

    return df


def plot_isoform_expression(
    pdf, 
    genes, 
    count_threshold=1000, 
    isoform_plot_kwargs={},
    gene_plot_kwargs={}
):
    """
    Plots isoform expression percentages as a stacked bar chart,
    with a secondary bar plot showing total gene expression (log CPM).

    Args:
        pdf (pandas.DataFrame): The dataframe containing expression data.
        genes (list): List of gene names to filter the analysis.
        count_threshold (int, optional): Minimum gene count to include. Defaults to 1000.
        isoform_plot_kwargs (dict, optional): Keyword arguments for the isoform bar plot.
        gene_plot_kwargs (dict, optional): Keyword arguments for the gene expression bar plot.
    """

    # Filter data
    table, pdf = pivot_isoforms(pdf, genes=genes, count_threshold=count_threshold)
    
    # Default plot kwargs if not provided
    default_isoform_kwargs = {
        'kind': 'bar',
        'stacked': True,
        'width': 0.75,
        'cmap': 'tab20b',
        'ec': 'k'
    }
    default_gene_kwargs = {
        'x': 'gene_name',
        'y': 'log_CPM',
        'ec': 'k',
        'width': 0.75,
        'color': 'lightgrey',
    }

    # Merge with provided kwargs (user-provided values override defaults)
    isoform_plot_kwargs = {**default_isoform_kwargs, **isoform_plot_kwargs}
    gene_plot_kwargs = {**default_gene_kwargs, **gene_plot_kwargs}

    # Plot isoform expression percentages
    ax1 = table.plot(**isoform_plot_kwargs)
    plt.ylabel("Percent of Expression")

    # Get total expression for the second plot
    pdf_gene_cpm = pdf[['gene_name', 'log_CPM']].drop_duplicates().sort_values(by='gene_name')

    # Add the second plot above the first
    ax_divider = make_axes_locatable(ax1)
    ax2 = ax_divider.append_axes("top", size="35%", pad="10%") 

    gene_plot_kwargs['ax'] = ax2
    sns.barplot(data=pdf_gene_cpm, **gene_plot_kwargs)

    ax2.sharex(ax1)
    ax2.tick_params(axis='x', labelsize=1, labelcolor='w')
    ax2.set_ylabel('CPM (log)')
    ax1.set_xlabel("")
    ax2.set_xlabel("")
    
    # Create a new axis for the legend
    ax_legend = ax_divider.append_axes("right", size="15%", pad="5%")
    ax_legend.axis("off")  # Turn off axis lines and labels for the legend

    # Move the legend to the new axis
    handles, labels = ax1.get_legend_handles_labels()
    ax_legend.legend(handles, 
                     labels, 
                     loc='center left',
                     title='Isoform')
    ax1.legend().remove()