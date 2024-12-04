# ONT-Single-Cell

This repository contains tools and workflows for processing and analyzing single-cell data generated with Oxford Nanopore Technologies (ONT). It is designed to support high-throughput analysis pipelines and flexible exploration of single-cell transcriptomic data.

## Features
- **Snakemake Workflow**: Automates data processing using a customizable pipeline defined in `Snakefile`.
- **Conda Environment**: Pre-configured environment setup for reproducible analyses.
- **Scripts and Utilities**: Includes helper scripts for specific tasks like data formatting, cluster configuration, and runtime management.
- **Resource Management**: SLURM-based launcher script for efficient use of computational resources.

## Repository Structure
- **`config/`**: Configuration files for Snakemake workflows.
- **`envs/`**: Conda environment specifications.
- **`notebooks/`**: Jupyter Notebooks for exploratory data analysis.
- **`resources/`**: Additional resources and reference data.
- **`rules/`**: Snakemake rule definitions.
- **`scripts/`**: Custom scripts used throughout the pipeline.
- **`utils/`**: Utility functions and helper tools.
- **`workflow.smk`**: The main workflow file for Snakemake.
