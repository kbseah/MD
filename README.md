# Metagenomic diversity index calculation

Snakemake workflow to calculate diversity indices from metagenomic predicted
proteome sequences.


## Usage

 * To run the pipeline, a recent version of `mamba` must be in the current
   path.
 * Specify sample aliases and paths to sequence files (Fasta format, a.a.
   sequences) in the config file `config/config.yaml`.
 * Run the script `run_pipeline.sh`: it will first check for a Conda
   environment containing Snakemake in the current folder at path
   `./snakemake-env`, otherwise it will create the environment and activate it
   before running the workflow.
 * Output will be written to `results/` and log fies to `logs/`


## Citation

[pending]
