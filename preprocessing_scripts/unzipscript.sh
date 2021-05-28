#!/bin/bash
FILES=$(find /data/caitken/eIF3sequencing/ -name "*.gz")

OIFS="$IFS"
IFS=$'\n'

for file in $FILES
do
#   NAME=$(basename $file .gz)
	NAME=$(echo $file | cut -d '/' -f 7 | cut -d '.' -f 1)
    gunzip  -kc $file > /data/caitken/files/${NAME}.fastq & 
    echo unzipping ${NAME}...
    wait
done

