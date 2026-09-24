#!/bin/bash
#SBATCH --job-name=HMSC_big_model
#SBATCH --account=nn9725k
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --array=0-4

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

input_path="init_GPP.preds.rds"
output_path=$(printf "post_chain%.2d_file.rds" $SLURM_ARRAY_TASK_ID)


apptainer exec --nv \
	hmsc-hpc-aarch64_0.1.8.sif \
	python3 -m hmsc.run_gibbs_sampler \
	--input "$input_path" \
	--output "$output_path"\
	--samples 500 \
	--thin 10 \
	--transient 20000 \
	--verbose 1000 \
	--chain "$SLURM_ARRAY_TASK_ID"


