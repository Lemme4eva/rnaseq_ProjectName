# RNA-Seq Analysis Project

This repository contains a reproducible RNA-Seq analysis workflow.

## Directory structure
- `data/raw_data/` — Raw sequencing files (FASTQ, etc.)
- `data/reference_data/` — Reference genome, annotation files
- `data/meta_data/` — Sample metadata
- `data/processed_data/` — Processed data: trimming, alignment, counts
- `results/` — QC, differential expression, and functional profiling results
- `reports/` — Analysis reports
- `scripts/` — R scripts for analysis steps
- `R/` — Internal functions
- `logs/` — Log files from pipelines

## Environment
This project uses **renv** for reproducibility.
Run `renv::restore()` to install the correct package versions.

## Version control
Git and GitHub are configured automatically by this setup script.
