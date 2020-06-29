#!/bin/bash

#SBATCH -J dZFW_cox
#SBATCH -o ./dZFW_cox.res
#SBATCH -e ./dZFW_cox.err

#SBATCH --mail-user laura.iacovissi@gmail.com
#SBATCH --mail-type=ALL

#SBATCH -c 15
#SBATCH --mem=10G
#SBATCH -p medium

# Run the python script
python3 DZFW_cox.py
