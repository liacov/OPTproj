#!/bin/bash

#SBATCH -J dZFW_lasso
#SBATCH -o ./dZFW_lasso.res
#SBATCH -e ./dZFW_lasso.err

#SBATCH --mail-user laura.iacovissi@gmail.com
#SBATCH --mail-type=ALL

#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH -p medium

# Run the python script
python3 DZFW_lasso.py
