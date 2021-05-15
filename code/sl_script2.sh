#!/bin/bash
#SBATCH -c 1     # Number of cores
#SBATCH --job-name=process_generated_data
#SBATCH -o process.out	# File to which STDOUT will be written
#SBATCH -e process.err	# File to which STDERR will be written
#SBATCH --output=process.txt
#SBATCH -p general		# Partition to submit to
#
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mail-user=astanciu@vassar.edu 
#SBATCH --mail-type=begin   # email me when the job starts 
#SBATCH --mail-type=end     # email me when the job finishes
#SBATCH --begin=05:00:00 # start at 5AM to ensure that the generator finished running.

squeue

python3 metagene_hypothesis.py --iterations 1000 simulated second_dataset.pickle simulated_data_ks_test_2.csv
