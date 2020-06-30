#!/bin/bash

#SBATCH -J sZFW_IRDSA_lasso
#SBATCH -o ./sZFW_IRDSA_lasso.res
#SBATCH -e ./sZFW_IRDSA_lasso.err

#SBATCH --mail-user laura.iacovissi@gmail.com
#SBATCH --mail-type=ALL

#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH -p medium

# Run the python script
python3 SZFW_lasso.py
