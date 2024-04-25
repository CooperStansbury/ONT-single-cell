import pandas as pd
import os
import sys


def get_input_names(inputs, output_dir):
    """A function to get output names from the basecall file """
    output_files = []
    for idx, row in inputs.iterrows():

        # handle differnet compression
        if row['file_path'].endswith('.gz'):
            ext = ".gz"
        else:
            ext = ''
            
        sample_id = row['sample_id']
        new_path = f"{output_dir}fastq/{sample_id}.raw.fastq{ext}"
        output_files.append(new_path)
    return output_files

