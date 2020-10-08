import pysam
import matplotlib.pyplot as plt
from pybedtools import BedTool
import pandas as pd
import sys

# This script is in ../code. Change these paths if needed
_kDefaultBamFile = "../AitkenLab/rnaseq_all.bam"
_kDefaultBedFile = "../indexes/yeast-all.bed"
_kDefaultMTIFFile = "../SteinmetzGilbert/stmtz_n_mtifs.txt"
_kDefaultChrom1Size = 230220
_kMaxReadLength = 100  # jic, but for chrI max read length is 44.

#gene:
#    - genome (sequence)
#    - genomic definition (interval - standard gene annotated and used in most RNAseq experiments)
#    - introns / exons (intervals)
#    - data on lengths (CDS, 5orf, 3urf)
#    - steinmetz data:
#        - Number of mTIF
#        - List of mTIFs
#        - Stats on tifs (sd, mean)

# TODO(astanciu) Download mRNA seq data from steinmetz

# Reads that violate certain isoforms.
# Maybe reduce the dataset to the reads that span the junction.


def grab_reads(chrom, reads, positive_strand_only=False):
    """
    Extracts all reads in chrom

    Parameters:
        chrom (string): name of chromosome
        reads (pysam.AlignmentFile): Experiment data

    Returns:
        list: List of gene objects 
    """

    filtered = reads.fetch(chrom)
    if positive_strand_only:
        return [x for x in filtered if not x.is_reverse]
    return list(reads.fetch(chrom))


def read_files(bedfile=_kDefaultBedFile,
               bamfile=_kDefaultBamFile,
               mtiffile=_kDefaultMTIFFile,
               positive_strand_only=False):
    bam = pysam.AlignmentFile(bamfile)
    # This is just an iterator in the file
    genes = BedTool(bedfile)  # Full genome loaded in memory
    mtifs = pd.read_csv(mtiffile,
                        sep='\t')  # Full pandas dataset loaded in memory

    if positive_strand_only:
        mtifs = mtifs[mtifs.strand == "+"]
    return bam, genes, mtifs


def generate_read_density(start, end, allreads):
    density = [0] * (end - start + 2)
    reads = grab_reads("chrI", allreads, positive_strand_only=True)
    for read in reads:
        # this condition checks whether a read overlaps at all with the interesting region [start, end]
        if read.reference_start <= end and read.reference_end >= start:
            read_adjusted_start = max(read.reference_start, start) - start
            read_adjusted_end = min(read.reference_end, end) - start
            density[read_adjusted_start] = density[read_adjusted_start] + 1
            density[read_adjusted_end + 1] = density[read_adjusted_end + 1] - 1
    for i in range(1, len(density)):
        density[i] = density[i] + density[i - 1]
    return density
