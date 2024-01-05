import pandas as pd
import numpy as np
import os
import sys

import pysam

        


if __name__ == "__main__":
    in_bam = sys.argv[1]
    out_bam = sys.argv[2]

    with pysam.AlignmentFile(in_bam, "rb") as bam_in:
        with pysam.AlignmentFile(out_bam, "wb", template=bam_in) as bam_out:
            for align in bam_in:
                barcode = align.query_name.split('_')[0]
                umi = align.query_name.split('_')[1].split("#")[0]
                read_name = align.query_name.split("#")[1][:-2]
    
                align.set_tag('CB', barcode, value_type="Z")
                align.set_tag('UB', umi, value_type="Z")
                align.set_tag('RD', read_name, value_type="Z")
    
    
                bam_out.write(align)


    