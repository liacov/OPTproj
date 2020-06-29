#!/bin/bash

#SBATCH -J inexactZFW_cox
#SBATCH -o ./inexactZFW_cox.res
#SBATCH -e ./inexactZFW_cox.err

#SBATCH --mail-user laura.iacovissi@gmail.com
#SBATCH --mail-type=ALL

#SBATCH -c 15
#SBATCH --mem=10G
#SBATCH -p medium

# Run the python script
python3 IZFW_cox.py
