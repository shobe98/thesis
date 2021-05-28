#!/bin/bash

OIFS="$IFS"
IFS=$'\n'

FILES=$(find /data/caitken/Ingolia/perfect -name "*.bam")

for f in $FILES
do
    echo "Indexing $f..."
    samtools index $f &
done


