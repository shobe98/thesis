#!/bin/bash
#SBATCH -c 1     # Number of cores
#SBATCH --job-name=data_generation
#SBATCH -o python.out	# File to which STDOUT will be written
#SBATCH -e python.err	# File to which STDERR will be written
#SBATCH --output=python_res.txt
#SBATCH -p general		# Partition to submit to
#
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mail-user=astanciu@vassar.edu 
#SBATCH --mail-type=begin   # email me when the job starts 
#SBATCH --mail-type=end     # email me when the job finishes


python3 generator.py
