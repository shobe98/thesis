#!/bin/bash

NAME=$(basename $1 .fq)

ROOTPATH="/data/caitken/Ingolia/intermediary/Riboseq"
idxpath="/data/caitken/indices"
FastQCpath="/data/caitken/Ingolia/FastQC_output"

mkdir -v -p $ROOTPATH 
mkdir -v -p $ROOTPATH

echo "Processing file $NAME"

#from original file to _clipped
fastx_clipper -a CTGTAGGCACCATCAAT -l 25 -c -n -v -i /data/caitken/files/${NAME}.fq -Q 33 -o ${ROOTPATH}/${NAME}_clipped.fq &
echo "removing adaptor..."
wait


#trimms the first nucleotide
#from _clipped to _trimmed
fastx_trimmer -Q 33 -i ${ROOTPATH}/${NAME}_clipped.fq -f 2 -o ${ROOTPATH}/${NAME}_trimmed.fq &
echo "trimming sequence..."
wait


# from _trimmed to _contaminating_rRNA and _rRNA_free
bowtie2 --seedlen=23 --un=${ROOTPATH}/${NAME}_rRNA_free.fq -x ${idxpath}/yeast_rRNA ${ROOTPATH}/${NAME}_trimmed.fq > ${ROOTPATH}/${NAME}_contaminating_rRNA.SAM &
echo "bowtie for removing rRNAs in riboseq reads"
wait

# starting tophat alignments

echo "starting tophat"

#from _rRNA_free to folder $NAME_aligned
#the folder contains the following:
#align_summary.txt
#accepted_hits.bam < the one we are interested in
#unmapped_hits.bam
#deletions.bed
#insertions.bed
#junctions.bed
# and log files
tophat --no-novel-juncs --GTF ${idxpath}/sac3.gtf -o ${ROOTPATH}/${NAME}_aligned ${idxpath}/saccharomyces_cerevisiae ${ROOTPATH}/${NAME}_rRNA_free.fq &
echo "generating alignment for ribos..."
echo "wait for 30 min or more"
wait

echo "file $NAME sucessfully processed"

echo "generating quality checks..."

mkdir -v -p ${FastQCpath}/logs
FQFILES=$(find $ROOTPATH -name "${NAME}*.fq")
for f in $FQFILES
do
    ~/FastQC/fastqc -o $FastQCpath > ${FastQCpath}/logs/${NAME}.out.log 2>  ${FastQCpath}/logs/${NAME}.err.log 
done

echo "removing unnecessary files"
ls -lR ${ROOTPATH}/${NAME}*
rm -f ${ROOTPATH}/${NAME}*_clipped.fq ${ROOTPATH}/${NAME}*_trimmed.fq ${ROOTPATH}/${NAME}*.SAM
mkdir -v -p ${FastQCpath}/logs
