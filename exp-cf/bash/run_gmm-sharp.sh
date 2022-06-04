#!/bin/sh
#
#SBATCH -A sml
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --nodelist=statler,waldorf
#SBATCH --mail-type=END

#call this from src
echo "Gamma is: ${SA}"

Rscript synthetic_exp-sharp.R --gmm_star ${SA} --alpha 0.2 --dtype ${DTYPE} --trial ${TRIAL} 

