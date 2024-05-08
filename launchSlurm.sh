#!/bin/bash

#SBATCH --mail-user=cstansbu@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=100G
#SBATCH --time=96:00:00
#SBATCH --account=indikar1
#SBATCH --partition=standard

SLURM_CONFIG='config/slurm'
CONFIG='config/config.yaml'
OUTPUT_PATH="$(cat ${CONFIG} | shyaml get-value output_path)"
THREADS=36

## build the workflow from the most current snakefile
cp Snakefile workflow.smk
echo "Built Workflow..."

snakemake --profile ${SLURM_CONFIG} --cores ${THREADS} --latency-wait 90 -s workflow.smk --use-conda --rerun-incomplete
