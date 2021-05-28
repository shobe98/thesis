#!/bin/bash

OIFS="$IFS"
IFS=$'\n'

FILES=$(find /data/caitken/Ingolia/intermediary/ -name "*.bam")

for f in $FILES
do
    FOLDER=$(dirname "$f")
    cd $FOLDER #this is important since we want our .bai files in the same folder as the .bam file
    pwd #print working directory (for extra safety)
    
    NAME=$(echo $f | cut -d '.' -f 1 | cut -d '/' -f 7)
    NAME="${NAME%_*}"
    
    OLD_NAME=$(basename $f)
    TYPE=$(basename $f .bam)
    TYPE="${TYPE##*_}"
    NEW_NAME="${NAME}_${TYPE}.bam"
    
    # rename all files from /lane_name/*.bam to
    # lane_name/lane_name_accepted.bam
    # or lane_name/lane_name_unmapped.bam
    # comment out if files are already properly named
    
    COMMAND="mv $OLD_NAME $NEW_NAME"    
    echo "command is: $COMMAND"
    eval "$COMMAND"
done

FILES=$(find /data/caitken/Ingolia/intermediary/ -name "*_accepted.bam")

for f in $FILES
do
    echo processing $f ...
    NAME="$(basename $f)"
    BASE="$(echo $NAME |  cut -d '.' -f 1)"
    samtools view -h $f | \
    grep -E '(NM:i:0)|(^@)' | samtools view -S -b - \
    > /data/caitken/Ingolia/perfect/${BASE}_perfect.bam
done



