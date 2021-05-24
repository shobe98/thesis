#!/bin/bash
#SBATCH -c 1     # Number of cores
#SBATCH --job-name=densities
#SBATCH -o densities.out	# File to which STDOUT will be written
#SBATCH -e densities.err	# File to which STDERR will be written
#SBATCH --output=densities.txt
#SBATCH -p general		# Partition to submit to
#
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mail-user=astanciu@vassar.edu 
#SBATCH --mail-type=begin   # email me when the job starts 
#SBATCH --mail-type=end     # email me when the job finishes


python3 generate_densities.py
