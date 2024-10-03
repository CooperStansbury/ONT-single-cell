import pandas as pd
import numpy as np
import os
import sys
import pysam

        
if __name__ == "__main__":
    in_bam = sys.argv[1]
    out_bam = sys.argv[2]
    out_record = sys.argv[3]

    records = []

    with pysam.AlignmentFile(in_bam, "rb") as bam_in:
        with pysam.AlignmentFile(out_bam, "wb", template=bam_in) as bam_out:
            for align in bam_in:
                barcode = align.query_name.split('_')[0]

                align.set_tag('NH', 1)
                bam_out.write(align)

                # structure the record table
                row = {
                    'barcode' : barcode,
                }

                records.append(row)

    # save the record table
    records = pd.DataFrame(records)
    records = records.drop_duplicates()
    records.to_csv(out_record, index=False, header=False)


    