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
 * Additional Snakemake parameters can be appended to the `run_pipeline.sh`
   command, e.g. `bash run_pipeline.sh --dryrun` for dry run (list workflow
   steps to be executed without actually running them).
 * Output will be written to `results/` and log files to `logs/`


## Citation

[pending]
