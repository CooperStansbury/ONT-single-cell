import pandas as pd
import numpy as np
import os
import sys
import pysam
import pyranges as pr


def load_gtf(gtf_path: str) -> pd.DataFrame:
    """
    Loads a GTF file and extracts relevant gene information for a specified list of transcription factors (TFs).

    This function reads a GTF file using the pr.read_gtf() function, filters for gene records, and further narrows 
    the selection to the specified TFs. It then aggregates information by gene name and feature, calculating the 
    minimum start, maximum end, and length for each gene.

    Args:
        gtf_path (str): The path to the GTF file.

    Returns:
        pd.DataFrame: A DataFrame containing the processed gene information, including columns for gene name, 
            chromosome, strand, start, end, and length.

    Example:
        genes_df = load_gtf("/path/to/your/gtf_file.gtf")
    """
    # Read the GTF file into a DataFrame
    gf = pr.read_gtf(gtf_path)
    gdf = gf.df

    # Define the transcription factors of interest
    TFs = ['GATA2', 'GFI1B', 'FOS', 'STAT5A', 'REL']

    # Filter for gene records and the specified TFs
    genes = gdf[gdf['Feature'].isin(['gene']) & gdf['gene_name'].isin(TFs)].copy()

    # Aggregate gene information, calculate length, and reset the index
    genes = genes.groupby(['gene_name', 'Feature']).agg(
        Chromosome=('Chromosome', 'first'),
        Strand=('Strand', 'first'),
        Start=('Start', 'min'),
        End=('End', 'max'),
    ).reset_index(drop=False)

    genes['Length'] = genes['End'] - genes['Start']

    return genes
        
    
def extract_query_names(bamfile: pysam.AlignmentFile, genes: pd.DataFrame, buffer_bp: int = 1000) -> pd.DataFrame:
    """
    Extracts query names (read IDs) from an already opened BAM file that overlap with specified gene regions,
    along with the corresponding gene information.

    This function iterates through gene records in a DataFrame, fetches reads from the BAM file
    that overlap with the gene regions (including a buffer), and extracts the query names of
    these reads. It then returns a DataFrame containing the query names and the corresponding
    gene information (gene name).

    Args:
        bamfile (pysam.AlignmentFile): An opened pysam.AlignmentFile object representing the BAM file.
        genes (pd.DataFrame): DataFrame containing gene information with columns 'Chromosome', 'Start', 'End', and 'gene_name'.
        buffer_bp (int, optional): Buffer size in base pairs to extend gene regions for fetching reads. Defaults to 1000.

    Returns:
        pd.DataFrame: A DataFrame containing columns 'query_name', 'gene_name'.
    """

    qnames = []

    for _, gene_rec in genes.iterrows():
        chrom = gene_rec['Chromosome']
        start = gene_rec['Start'] - buffer_bp
        end = gene_rec['End'] + buffer_bp

        for read in bamfile.fetch(chrom, start, end):
            qnames.append({
                'query_name': read.query_name,
                'gene_name': gene_rec['gene_name']
            })

    return pd.DataFrame(qnames)


if __name__ == "__main__":
    in_bam = sys.argv[1]
    gtf_path = sys.argv[2]
    read_id_outfile = sys.argv[3]
    read_map_outfile = sys.argv[4]

    # load the reprogramming genes
    genes = load_gtf(gtf_path)

    buffer_bp = 1000 # base pair fudge factor    
    bamfile = pysam.AlignmentFile(in_bam, "rb")
    
    qnames = extract_query_names(bamfile, genes, buffer_bp)
    
    # export unique read names
    unique_queries = qnames['query_name'].unique()
    read_ids = pd.DataFrame({'qname' : unique_queries})
    read_ids.to_csv(read_id_outfile, index=False, header=False,)
    
    # export 
    qnames.to_csv(read_map_outfile,index=False,)
    
    

    

    


    