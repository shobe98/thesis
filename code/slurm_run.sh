#!/bin/bash
#SBATCH -c 1     # Number of cores
#SBATCH --job-name=thesis_run_all_chrs
#SBATCH -o my_job_%j.out	# File to which STDOUT will be written
#SBATCH -e errors_%j.err	# File to which STDERR will be written
#SBATCH --output=res.txt
#SBATCH -p general		# Partition to submit to
#
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL  	  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=astanciu@vassar.edu

export MPLBACKEND=Agg

yes | pip3 install pysam
yes | pip3 install matplotlib
yes | pip3 install pandas
yes | pip3 install pybedtools
yes | pip3 install numpy

python3 tester.py
