#!/bin/sh
#SBATCH --mem=200g
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J geneLevelOnly
#SBATCH -a 1-10
#SBATCH -o /data/smithral/sgacc_wgcna/RCCA_toolkit/cca_pls_toolkit_final/logs/log_%A_%a.out
#SBATCH --mail-type=ALL

module purge
module load matlab

mu=(0.1)

matlab -nodisplay -nodesktop -nojvm -nosplash<<EOF
	cd /data/smithral/sgacc_wgcna/RCCA_toolkit/cca_pls_toolkit_dev_final
	grcca_analysis($mu)	
	exit
EOF
