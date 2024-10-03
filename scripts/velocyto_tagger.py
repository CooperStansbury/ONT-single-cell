import pandas as pd
import numpy as np
import os
import sys
import time
import pysam

        
if __name__ == "__main__":
    in_bam = sys.argv[1]
    out_bam = sys.argv[2]

    print(f"Input BAM file: {in_bam}")
    print(f"Output BAM file: {out_bam}")

    records = []  # This list seems unused, consider removing it

    start_time = time.time()
    processed_reads = 0
    update_interval = 100000  # Adjust this for how often you want updates

    with pysam.AlignmentFile(in_bam, "rb") as bam_in:
        with pysam.AlignmentFile(out_bam, "wb", template=bam_in) as bam_out:
            for align in bam_in:
                barcode = align.query_name.split('_')[0]
                umi = align.query_name.split('_')[1].split("#")[0]
                read_name = align.query_name.split("#")[1][:-2]

                align.set_tag('CB', barcode, value_type="Z")
                align.set_tag('UB', umi, value_type="Z")
                align.set_tag('RD', read_name, value_type="Z")
                align.set_tag('NH', 1)
                bam_out.write(align)

                processed_reads += 1
                if processed_reads % update_interval == 0:
                    elapsed_time = time.time() - start_time
                    print(f"Processed {processed_reads} reads "
                          f"({elapsed_time:.2f} seconds elapsed)")

    end_time = time.time()
    total_time = end_time - start_time

    print("Done!")
    print(f"Total reads processed: {processed_reads}")
    print(f"Total time taken: {total_time:.2f} seconds")

            


    