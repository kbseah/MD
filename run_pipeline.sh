#!/bin/bash

eval "$(conda shell.bash hook)"
# Test for snakemake environment in current folder
mamba list -p ./snakemake-env &> /dev/null
RETVAL=$?
if [[ $RETVAL == 1 ]]
then
  echo "Creating environment with Snakemake at ./snakemake-env"
  mamba env create -p ./snakemake-env --file "$PWD/workflow/envs/snakemake.yml"
else
  echo "Environment with Snakemake found at ./snakemake-env"
fi
CORES=12
echo "Activating Conda environment"
conda activate ./snakemake-env
snakemake --sdm conda --conda-frontend mamba --cores ${CORES} \
  --configfile config/config.yaml $@
