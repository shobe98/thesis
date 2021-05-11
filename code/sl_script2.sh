#!/bin/bash
#SBATCH -c 1     # Number of cores
#SBATCH --job-name=data_generation
#SBATCH -o simulated.out	# File to which STDOUT will be written
#SBATCH -e simulated.err	# File to which STDERR will be written
#SBATCH --output=simulated_res.txt
#SBATCH -p general		# Partition to submit to
#
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mail-user=astanciu@vassar.edu 
#SBATCH --mail-type=begin   # email me when the job starts 
#SBATCH --mail-type=end     # email me when the job finishes

python3 metagene_hypothesis.py --iterations 1000 simulated first_dataset.pickle simulated_data_ks_test.csv
