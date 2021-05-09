#!/bin/bash
#SBATCH -c 1     # Number of cores
#SBATCH --job-name=bad_genes
#SBATCH -o logs/bad_genes_%j.out	# File to which STDOUT will be written
#SBATCH -e logs/bad_genes_%j.err	# File to which STDERR will be written
#SBATCH --output=logs/bad_genes_%j.txt
#SBATCH -p general		# Partition to submit to
#
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mail-user=astanciu@vassar.edu 
#SBATCH --mail-type=begin   # email me when the job starts 
#SBATCH --mail-type=end     # email me when the job finishes


python3 code/identify_genes_with_too_many_reads.py
