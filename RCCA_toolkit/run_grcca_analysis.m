#!/bin/bash
#SBATCH --job-name=grcca_mu_array
#SBATCH --output=/data/smithral/sgacc_wgcna/RCCA_toolkit/cca_pls_toolkit_final/logs/log_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200g
#SBATCH --array=1-10
#SBATCH --mail-type=ALL

# Clean environment and load matlab
module purge
module load matlab

# Define the range of mu values using awk for logarithmic spacing
mu_values=($(awk 'BEGIN{for(i=0;i<=9;i++) printf "%.4f ", 10^(log(0.1)/log(10) + i*(log(10)/log(10)/9))}'))

# Get the mu value for the current job
mu=${mu_values[$SLURM_ARRAY_TASK_ID-1]}

# Run MATLAB with the current mu value
matlab -nodisplay -nodesktop -nojvm -nosplash<<EOF
    cd /data/smithral/sgacc_wgcna/RCCA_toolkit/cca_pls_toolkit_dev_final
    grcca_analysis($mu)
    exit
EOF

# Print a message indicating completion
echo "Completed analysis for mu = $mu"