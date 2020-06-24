#!/bin/bash

#SBATCH -J sZFW_IRDSA_cox
#SBATCH -o ./sZFW_IRDSA_cox.res
#SBATCH -e ./sZFW_IRDSA_cox.err

#SBATCH --mail-user laura.iacovissi@gmail.com
#SBATCH --mail-type=ALL

#SBATCH --cpus-per-task=24
#SBATCH --mem=10G
#SBATCH --partition=short

# Run the python script
python3 SZFW_cox.py
