import pandas as pd
import numpy as np
import os
import sys

def create_read_map(read_map_path: str) -> dict[str, str]:
    """Creates a dictionary mapping query names to gene names from a CSV file.

    Args:
        read_map_path (str): The file path to the CSV containing the mapping.

    Returns:
        dict: A dictionary where keys are query names and values are gene names.
    """
    
    read_map_df = pd.read_csv(read_map_path)
    return read_map_df.set_index('query_name')['gene_name'].to_dict()

def parse_seqID(input_string):
    """
    Parse an input string containing barcode, UMI, and read name.

    Args:
        input_string (str): A string containing barcode, UMI, and read name separated by '#' and '_'.

    Returns:
        tuple: A tuple containing barcode, UMI, and read name.
    """
    parts = input_string.split('#')
    barcode_umi = parts[0]
    barcode, umi = barcode_umi.split('_')
    read_name = parts[1].split('_')[0]
    return barcode, umi, read_name


def process_read_metadata(tmp: pd.DataFrame) -> pd.DataFrame:
    """Processes read metadata to identify V5 and HT patterns in gene groups.

    Args:
        tmp (pd.DataFrame): Input DataFrame containing read metadata, 
                            including columns 'read_name', 'barcode', 'umi', 'strand', 'gene_name', and 'patternName'.

    Returns:
        pd.DataFrame: A DataFrame summarizing the presence of V5 and HT patterns for each unique combination of 
                      read name, barcode, UMI, strand, and gene name.

    """
    res = []

    for read_meta, group in tmp.groupby(['read_name', 'barcode', 'umi', 'strand', 'gene_name']):
        read_name, barcode, umi, strand, gene_name = read_meta

        has_V5 = 'V5' in group['patternName'].str.strip().values
        has_HT = 'HT' in group['patternName'].str.strip().values

        row = {
            'barcode': barcode,
            'umi': umi,
            'read_name': read_name,
            'gene_name': gene_name,
            'has_V5': has_V5,
            'has_HT': has_HT,
        }
        res.append(row)

    return pd.DataFrame(res)


def count_has_vh_by_unique_umi(df):
    """Counts the number of has_V5 and has_HT columns by unique UMIs per barcode per gene_name.

    Args:
    df: A pandas DataFrame with columns 'barcode', 'umi', 'read_name', 'gene_name',
        'has_V5', and 'has_HT'.

    Returns:
    A pandas DataFrame with columns 'barcode', 'gene_name', 'num_has_V5', 'num_has_HT',
        where 'num_has_V5' is the count of True values in the 'has_V5' column for
        each unique combination of 'barcode' and 'gene_name', and 'num_has_HT' is
        similarly defined for the 'has_HT' column.
    """

    return (df.groupby(['barcode', 'gene_name'])
          .agg(num_V5=('has_V5', sum), num_HT=('has_HT', sum))).reset_index()


def process_barcode_search(file_list, read_map):
    """
    Processes sequencing files to extract barcode information, map reads to genes, and
    summarize the presence of V5 and HT tags.

    This function iterates through a list of sequencing files, parses the sequence IDs to
    extract barcode, UMI, and read name information. It then maps reads to their corresponding
    gene names and processes read metadata to identify V5 and HT tags. Finally, it aggregates
    the results into a DataFrame.

    Args:
        file_list (list): A list of file paths to sequencing data files.
        read_map (dict): A dictionary mapping read IDs to their corresponding gene names.

    Returns:
        pandas.DataFrame: A DataFrame containing the processed data, including columns for
            barcode, gene name, and counts of V5 and HT tags.

    Example:
        result_df = process_barcode_search(["file1.txt", "file2.txt"], read_map)
    """

    df = []

    for fpath in file_list:
        file_id = os.path.basename(fpath).split(".")[0]

        # load data 
        tmp = pd.read_csv(fpath, sep="\t")
        if tmp.shape[0] == 0:
            continue

        tmp[['barcode', 'umi', 'read_name']] = tmp['seqID'].apply(parse_seqID).apply(pd.Series)

        # map the reads to the genes that they aligned to
        tmp['gene_name'] = tmp['seqID'].map(read_map)

        # create boolean flags for each detected V5 tag
        tmp = process_read_metadata(tmp)
        tmp = count_has_vh_by_unique_umi(tmp)
        tmp['file_id'] = file_id
        df.append(tmp)

    df = pd.concat(df)
    return df


def calculate_sum_by_gene_barcode(df):
    """
    Calculates the sum of 'num_V5' and 'num_HT' for each gene by barcode,
    aggregating over all 'file_id'.

    Args:
        df: A pandas DataFrame with columns 'barcode', 'gene_name', 'num_V5', 'num_HT', and 'file_id'.

    Returns:
        A pandas DataFrame with columns 'barcode', 'gene_name', 'num_V5', and 'num_HT', 
        where 'num_V5' and 'num_HT' are the summed values.
    """

    # Group by 'barcode' and 'gene_name', then sum the relevant columns
    result_df = df.groupby(['barcode', 'gene_name'])[['num_V5', 'num_HT']].sum().reset_index()
    
    result_df = pd.pivot_table(result_df, index='barcode',
                               columns='gene_name',
                               fill_value=0)

    result_df.columns = ['-'.join(col).strip().replace("num_", "") for col in result_df.columns.values]
    return result_df


        
if __name__ == "__main__":
    read_map_path = sys.argv[1]
    out_path = sys.argv[2]
    file_list = sys.argv[3:]
    
    
    # create read_map
    read_map = create_read_map(read_map_path)
    
    # load all files and aggregate
    df = process_barcode_search(file_list, read_map)
    df = calculate_sum_by_gene_barcode(df)
    
    df = df.reset_index()
    df.to_csv(out_path, index=False)
    
    
    