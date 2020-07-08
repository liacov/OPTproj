#!/bin/bash

#SBATCH -J inexactZFW_lasso_mind
#SBATCH -o ./inexactZFW_lasso_mind.res
#SBATCH -e ./inexactZFW_lasso_mind.err

#SBATCH --mail-user laura.iacovissi@gmail.com
#SBATCH --mail-type=ALL

#SBATCH -c 10
#SBATCH --mem=10G
#SBATCH -p long

# Run the python script
python3 IZFW_lasso_long.py
