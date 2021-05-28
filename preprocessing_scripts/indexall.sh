#!/bin/bash

OIFS="$IFS"
IFS=$'\n'

FILES=$(find /data/caitken/ribo/ -name "*.bam")

for f in $FILES
do
    FOLDER=$(dirname "$f")
    cd $FOLDER #this is important since we want our .bai files in the same folder as the .bam file
    pwd #print working directory (for extra safety)
    
    # the file is in a path that looks like this /data/caitken/some/other/folders/the_name_of_the_processed_file_aligned/aligned.bam
    # or /data/caitken/some/other/folders/the_name_of_the_processed_file_aligned/unmapped.bam

    # take the name of the file, split it in two arounf '.'. Take the string before the extension, split it by '/' and take the 6th string. This should be changed and tested according to your own folder structure. 
    NAME=$(echo $f | cut -d '.' -f 1 | cut -d '/' -f 6)
    NAME="${NAME%_*}"
    
    OLD_NAME=$(basename $f)
    TYPE=$(basename $f .bam)
    TYPE="${TYPE%_*}"
    NEW_NAME="${NAME}_${TYPE}.bam"
    
    # rename all files from /lane_name/*.bam to
    # lane_name/lane_name_accepted.bam
    # or lane_name/lane_name_unmapped.bam
    # comment out if files are already properly named

    COMMAND="mv $OLD_NAME $NEW_NAME"    
    echo "command is: $COMMAND"
    eval "$COMMAND"

    COMMAND="samtools index $NEW_NAME"
    echo "command is: $COMMAND"
    eval "$COMMAND"

    echo "done indexing $f"
done
