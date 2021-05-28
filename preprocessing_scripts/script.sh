#!/bin/bash

echo "unzipping scripts... already done!"
#./unzipscript.sh
wait

echo "renaming files... already done!"
#./rename.sh
#wait

#module load bowtie #uncomment for running on the cluster
#module load tophat #uncomment for running on the cluster

ROOTPATH="/data/caitken/Ingolia"
idxpath="/data/caitken/indices"

mkdir -v -p $ROOTPATH


bowtie2-build ${idxpath}/yeast_rRNA.fa ${idxpath}/yeast_rRNA &
echo "indexing rRNA..."
wait

#indexing saccharomyces_cerevisiae with bowtie, file needs to be downloaded
bowtie2-build ${idxpath}/sacCer3.fa ${idxpath}/saccharomyces_cerevisiae &
echo "indexing scchromyces..."
wait

OIFS="$IFS"
IFS=$'\n'

RIBOS=$(find /data/caitken/files/ -name "*CA?_1*.fq"; find /data/caitken/files/ -name "*CA?_2*.fq")

for f in $RIBOS
do
    ./riboscript.sh $f
done

echo "finished alignments for ribo protected sequences and bam files are generated in _aligned folders"

MRNAS=$(find /data/caitken/files/ -name "*CA?_3*.fq"; find /data/caitken/files/ -name "*CA?_4*.fq")

for f in $MRNAS
do
    ./mrnascript.sh $f
done

echo "finished alignments for mrna sequences and bam files are generated in _aligned folders"

