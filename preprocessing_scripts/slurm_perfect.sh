#!/bin/bash
echo "Perfect reads"
./perfect_reads.sh > ../logs/perfect.out 2> ../logs/perfect.err 
echo "Running cksum"
cksum /data/caitken/Ingolia/perfect/* > /data/caitken/Ingolia/perfect/checksums.txt
echo "Done"
