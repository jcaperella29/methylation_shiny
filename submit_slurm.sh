#!/bin/bash
#SBATCH --job-name=shiny_methylation
#SBATCH --output=shiny_app.log
#SBATCH --error=shiny_app.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --partition=debug

# Run container
singularity exec methylation_shiny.sif Rscript /app/start.R
