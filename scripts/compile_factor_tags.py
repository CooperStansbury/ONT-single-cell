import pandas as pd
import numpy as np
import os
import sys


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


def read_and_combine_files(file_list):
    """
    Reads and combines multiple tab-separated files containing sequencing data.

    This function efficiently processes each file in the `file_list`, parsing the 'seqID' column 
    into 'barcode', 'umi', and 'read_name' using the `parse_seqID` function defined elsewhere in 
    the file. It then selects the relevant columns and concatenates the data into a single Pandas 
    DataFrame. Empty files are automatically skipped.

    Args:
        file_list (list): A list of file paths to read.

    Returns:
        pandas.DataFrame: A combined DataFrame containing the parsed and filtered data.
            The DataFrame will have columns: 'barcode', 'umi', 'read_name', 'patternName', 
            'strand', 'start', and 'end'. If all input files are empty, an empty DataFrame
            with the correct columns will be returned.
    """

    df = pd.concat(
        (
            pd.read_csv(file_path, sep="\t")
            .assign(
                barcode=lambda df: df["seqID"].apply(parse_seqID).str[0],
                umi=lambda df: df["seqID"].apply(parse_seqID).str[1],
                read_name=lambda df: df["seqID"].apply(parse_seqID).str[2],
            )
            [[
                "barcode",
                "umi",
                "read_name",
                "patternName",
                "strand",
                "start",
                "end",
            ]]
            for file_path in file_list
            if not pd.read_csv(file_path, sep="\t").empty
        )
    )
    return df

        
if __name__ == "__main__":
    out_path = sys.argv[1]
    file_list = sys.argv[2:]
    
    df = read_and_combine_files(file_list)
    df.to_csv(out_path, index=False)
    
    
    