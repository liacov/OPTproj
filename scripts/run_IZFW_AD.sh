#!/bin/bash

#SBATCH -J inexactZFW_AD_long
#SBATCH -o ./inexactZFW_AD_long.res
#SBATCH -e ./inexactZFW_AD_long.err

#SBATCH --mail-user laura.iacovissi@gmail.com
#SBATCH --mail-type=ALL

#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH -p long

# Run the python script
python3 IZFW_AD_long.py
