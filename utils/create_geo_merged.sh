#!/bin/bash

# Define input and output directories
input_dir="/scratch/indikar_root/indikar1/cstansbu/HSC/fastq"
output_dir="/nfs/turbo/umms-indikar/shared/projects/HSC/geo_submission/data"

# Loop through all .fastq.gz files in the input directory
for file in "$input_dir"/*.fastq.gz; do
  # Extract the flow cell ID from the first line of the file
  flow_cell_id=$(zcat "$file" | head -n 1 | sed 's/.*flow_cell_id=\([^ ]*\).*/\1/')

  # Check if the flow cell ID starts with "F"
  if [[ "$flow_cell_id" == F* ]]; then
    # Print source, flow cell ID, and destination
    printf "File: %s\nFlow cell ID: %s\nWriting to: %s/gridion_ihsc.fastq.gz\n\n" "$file" "$flow_cell_id" "$output_dir"
    # Write the file to gridion_ihsc.fastq.gz in the output directory
    zcat "$file" >> "$output_dir"/gridion_ihsc.fastq.gz
  else
    # Print source, flow cell ID, and destination
    printf "File: %s\nFlow cell ID: %s\nWriting to: %s/p2_ihsc.fastq.gz\n\n" "$file" "$flow_cell_id" "$output_dir"
    # Write the file to p2_ihsc.fastq.gz in the output directory
    zcat "$file" >> "$output_dir"/p2_ihsc.fastq.gz
  fi
done