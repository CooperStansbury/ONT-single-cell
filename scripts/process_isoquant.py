import pandas as pd
import numpy as np
import os
import sys
import pyranges as pr

def load_gene_table(fpath: str = "/scratch/indikar_root/indikar1/cstansbu/HSC/references/geneTable.csv") -> pd.DataFrame:
    """
    Loads a gene table CSV file from the specified path, 
    cleans the data, and returns a Pandas DataFrame.

    Args:
        fpath: The path to the CSV file.

    Returns:
        A cleaned pandas DataFrame containing the gene data.
    """

    columns = [
        'transcript_id',
        'transcript_name',
        'transcript_biotype',
        'gene_id',
        'gene_name',
        'gene_biotype',
    ]

    # Read only the necessary columns
    df = pd.read_csv(fpath, low_memory=False, usecols=columns)

    # Drop rows with NaN values in key columns and remove duplicates
    df = df.dropna(subset=["transcript_id", "gene_id"]).drop_duplicates()
    return df


def map_counts_to_df(df, gene_path, transcript_path):
    """Maps gene and transcript counts from CSV files to a DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame containing 'gene_id' and 'transcript_id' columns.
        gene_path (str): Path to the CSV file containing gene counts.
        transcript_path (str): Path to the CSV file containing transcript counts.

    Returns:
        pd.DataFrame: The input DataFrame with additional columns 'gene_counts' and 'transcript_counts'.
    """
    # Load gene counts
    genes_df = pd.read_csv(gene_path, usecols=['gene_id', 'total_counts'])
    gene_counts = dict(zip(genes_df['gene_id'], genes_df['total_counts']))

    # Load transcript counts
    transcript_df = pd.read_csv(transcript_path, usecols=['gene_id', 'total_counts'])
    transcript_counts = dict(zip(transcript_df['gene_id'], transcript_df['total_counts']))

    # Map counts to the DataFrame
    df['gene_counts'] = df['gene_id'].map(gene_counts)
    df['gene_counts'] = df['gene_counts'].astype(int)
    df['transcript_counts'] = df['transcript_id'].map(transcript_counts)
    df['transcript_counts'] = df['transcript_counts'].astype(int)

    return df


def calculate_gene_cpm(df):
    """Calculates CPM (Counts Per Million) for genes from a DataFrame.

    Args:
        df (pd.DataFrame): A DataFrame containing columns 'gene_name' and 'gene_counts'.

    Returns:
        dict: A dictionary mapping gene names to their calculated CPM values.
    """

    # Filter and remove duplicates
    t = df[['gene_name', 'gene_counts']].drop_duplicates()

    # Remove rows with missing gene names
    t = t[t['gene_name'].notna()]

    # Calculate CPM
    total_counts = t['gene_counts'].sum()
    t['CPM'] = 1e6 * (t['gene_counts'] / total_counts)

    # Create gene to CPM dictionary
    gene_cpm = dict(zip(t['gene_name'], t['CPM']))

    return gene_cpm


        
if __name__ == "__main__":
    gene_table_path = sys.argv[1]
    gene_path = sys.argv[2]
    transcript_path = sys.argv[3]
    outpath = sys.argv[4]
    
    # load the gene information 
    df = load_gene_table(gene_table_path)
    
    # map the counts
    df = map_counts_to_df(df, gene_path, transcript_path)
    
    # compute and merge cpm
    gene_tpm = calculate_gene_cpm(df)
    df['gene_CPM'] = df['gene_name'].map(gene_tpm)
    
    # save the output
    df.to_csv(outpath, index=False)
    
    


   