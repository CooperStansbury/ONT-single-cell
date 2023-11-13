import sys
import pandas as pd
import numpy as np
import pyranges as pr
    
    
if __name__ == "__main__":
    gtfPath = sys.argv[1]
    outpath = sys.argv[2]
    
    # open gtf and save as csv
    gr = pr.read_gtf(gtfPath)
    df = gr.as_df()
    
    df.to_csv(outpath, index=False)