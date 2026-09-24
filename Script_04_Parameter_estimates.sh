#!/bin/bash
#SBATCH --job-name=Parameter_estimates_HMSC_model
#SBATCH --account=nn9725k
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --partition=large
#SBATCH --ntasks=1


## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module purge

module load NRIS/CPU
module load R/4.5.2-gfbf-2025b


R < Parameter_estimates.R --no-save
