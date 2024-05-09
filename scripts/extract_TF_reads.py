import pandas as pd
import numpy as np
import os
import sys
import pysam
import pyranges as pr

def load_gtf(gtf_path):
    # get the gtf records we care about
    gf = pr.read_gtf(gtf_path)
    gdf = gf.df
    
    TFs = [
        'GATA2', 
        'GFI1B', 
        'FOS', 
        'STAT5A',
        'REL',  
    ]
    
    print(f"{len(TFs)=}")
    
    genes = gdf.copy()
    genes = gdf[gdf['Feature'].isin(['gene'])]
    genes = genes[genes['gene_name'].isin(TFs)]
    
    genes = genes.groupby(['gene_name', 'Feature']).agg(
        Chromosome = ('Chromosome', 'first'),
        Strand = ('Strand', 'first'),
        Start = ('Start', 'min'),
        End = ('End', 'max'),
    ).reset_index(drop=False)
    
    genes['Length'] = genes['End'] - genes['Start']
    return genes
        


if __name__ == "__main__":
    in_bam = sys.argv[1]
    gtf_path = sys.argv[2]
    out_file = sys.argv[3]

    # load the reprogramming genes
    genes = load_gtf(gtf_path)

    buffer_bp = 1000 # base pair fudge factor    
    bamfile = pysam.AlignmentFile(in_bam, "rb")
    
    qnames = []
    for _, gene_rec in genes.iterrows():
        chrom = gene_rec['Chromosome']
        start = gene_rec['Start'] - buffer_bp
        end = gene_rec['End'] + buffer_bp
        
        for read in bamfile.fetch(chrom, start, end):
            query_name = read.query_name
            row = {
                'query_name' : query_name,
            }
            qnames.append(row)
        
    qnames = pd.DataFrame(qnames)
    qnames.to_csv(out_file, index=False, header=False,)

    

    


    