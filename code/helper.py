import pysam
import matplotlib.pyplot as plt
from pybedtools import BedTool
import pandas as pd
import sys
from numpy import unique
import numpy as np

# This script is in ../code. Change these paths if needed
_kDefaultBamFile = "../AitkenLab/rnaseq_all.bam"
_kDefaultBedFile = "../indexes/yeast-all.bed"
_kDefaultMTIFFile = "../SteinmetzGilbert/stmtz_n_mtifs.txt"
_kDefaultTIFsFile = "../parsed_steinmetz_s1_tifs.txt"
_kDefaultChrom1Size = 230220
_kMaxReadLength = 100  # jic, but for chrI max read length is 44.

kYeastChroms = [
    "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
    "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV"
]

#gene:
#    - genome (sequence)
#    - genomic definition (interval - standard gene annotated and used in most RNAseq experiments)
#    - introns / exons (intervals)
#    - data on lengths (CDS, 5orf, 3urf)
#    - steinmetz data:
#        - Number of mTIF
#        - List of mTIFs
#        - Stats on tifs (sd, mean)


class Gene:
    def __init__(
            self,
            sequence,  # 
            introns,
            exons,
            gene_range,
            cds,  # use gene_range and cds to get utr lengths
            steinmetz):
        self.sequence = sequence
        self.introns = introns
        self.exons = exons
        self.gene_range = gene_range
        self.cds = cds
        self.steinmetz = steinmetz


# TODO(astanciu) Download mRNA seq data from steinmetz


# Reads that violate certain isoforms.
# Maybe reduce the dataset to the reads that span the junction.
def read_all_tifs(fname=_kDefaultTIFsFile,
                  positive_strand_only=False,
                  return_dict=False):
    tifs = pd.read_table(fname)
    # Only keep tifs that overlap one uorf
    tifs = tifs[tifs.type == "Covering one intact ORF"]
    if positive_strand_only:
        tifs = tifs[tifs.strand == "+"]

    if not return_dict:
        return tifs

    grouped_tifs = {}
    for chrom in kYeastChroms:
        grouped_tifs[chrom] = read_tifs_chrom(tifs, chrom)
    return grouped_tifs


def read_tifs_chrom(tifs, chrom, positive_strand_only=False):
    if isinstance(chrom, str):
        chrom = kYeastChroms.index(chrom) + 1
    chrom_tifs = tifs[tifs.chr == chrom]
    if positive_strand_only:
        chrom_tifs = chrom_tifs[chrom_tifs.strand == "+"]

    genes = list(unique(chrom_tifs.name))
    grouped_tifs = {}
    for gene in genes:
        grouped_tifs[gene] = chrom_tifs[chrom_tifs.name == gene]
    return grouped_tifs


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


def generate_5utr_isoform_starts(tifs):
    # This is a naive version for now
    return list(unique(tifs.t5))


def save_density_plot(density, ox, filename):
    if not ox:
        ox = range(0, len(density))
    plt.plot(ox, density)
    plt.savefig(filename)
    plt.close()


def generate_metagene(reads, genes, tifs, chrom):
    chrom_tifs = read_tifs_chrom(tifs, "chrI", positive_strand_only=True)
    density = generate_read_density(0, _kDefaultChrom1Size, reads, genes)
    #TODO(astanciu): implement a funciton that generates the density of reads for a whole chromosome (no star/end params)

    splits = []

    for gene, tifs in chrom_tifs.items():
        split_positions = generate_5utr_isoform_starts(tifs)
        splits = splits + split_positions
        #list(np.random.choice(split_positions, 1, replace=False))
    make_metagene_plot(density, splits, "meta_chr1.png", 200, 200)


def generate_random_metagene(reads, genes, chrom, sample_size):
    density = generate_read_density(0, _kDefaultChrom1Size, reads, genes)
    splits = np.random.randint(201,
                               _kDefaultChrom1Size - 201,
                               size=sample_size)
    make_metagene_plot(density, splits, "meta_chr1_random_sample.png", 200,
                       200)


def make_metagene_plot(density, splits, plotname, upstream, downstream):
    print("Number of split sites is: " + str(len(splits)))
    meta = np.zeros(upstream + downstream)
    for splitidx in splits:
        meta = np.add(meta,
                      density[(splitidx - upstream):(splitidx + downstream)])

    save_density_plot(meta, range(-200, 200), plotname)


def read_files(bedfile=_kDefaultBedFile,
               bamfile=_kDefaultBamFile,
               mtiffile=_kDefaultMTIFFile,
               positive_strand_only=False):
    bam = pysam.AlignmentFile(bamfile)
    # This is just an iterator in the file
    genes = BedTool(bedfile)  # Full genome loaded in memory
    genes_dict = {}
    for g in genes:
        genes_dict[g.name] = g

    mtifs = pd.read_csv(mtiffile,
                        sep='\t')  # Full pandas dataset loaded in memory

    if positive_strand_only:
        mtifs = mtifs[mtifs.strand == "+"]
    return bam, genes_dict, mtifs


def generate_read_density(start, end, allreads, genome):
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
