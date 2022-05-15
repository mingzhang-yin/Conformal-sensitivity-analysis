#!/bin/sh
#
#SBATCH -A sml
#SBATCH --cpus-per-task=2
#SBATCH -t 10:00:00
##SBATCH --mail-user=js5334@columbia.edu
##SBATCH --mail-type=END

#call this from src
echo "Gamma is: ${SA}"

Rscript exp2_fct_loop.R --gmm_star ${SA} --alpha 0.2 --dtype ${DTYPE}


