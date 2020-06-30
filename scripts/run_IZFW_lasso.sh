#!/bin/bash

#SBATCH -J inexactZFW_lasso
#SBATCH -o ./inexactZFW_lasso.res
#SBATCH -e ./inexactZFW_lasso.err

#SBATCH --mail-user laura.iacovissi@gmail.com
#SBATCH --mail-type=ALL

#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH -p medium

# Run the python script
python3 IZFW_lasso.py
