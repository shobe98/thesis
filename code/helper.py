import pysam
import matplotlib.pyplot as plt
from pybedtools import BedTool
import pandas as pd
import sys
from numpy import unique
import numpy as np

# 213 Days til graduation
np.random.seed(213)

# This script is in ../code. Change these paths if needed
_kDefaultBamFile = "../AitkenLab/rnaseq_all.bam"
_kDefaultBedFile = "../indexes/yeast-all.bed"
_kDefaultMTIFFile = "../SteinmetzGilbert/stmtz_n_mtifs.txt"
_kDefaultTIFsFile = "../parsed_steinmetz_s1_tifs.txt"

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
    """
    Reads all TIFs from steinmetz (parsed) dataset.

    Parameters:
        fname (string): File name
        positive_strand_only (bool): whether to filter out tifs on negative strand
        return_dict (bool): if false returns whole dataset, otherwise returns a nested dict {chrom: gene: tifs(pandas)}

    Returns:
        Either a pandas dataset of all tifs or a dict {chrom: gene: tifs(pandas)}
    """
    tifs = pd.read_table(fname)
    # Only keep tifs that overlap one uorf
    tifs = tifs[tifs.type == "Covering one intact ORF"]
    if positive_strand_only:
        tifs = tifs[tifs.strand == "+"]

    # TODO(astanciu): this seems like bad design change it into two separate functions
    if not return_dict:
        return tifs

    grouped_tifs = {}
    for chrom in kYeastChroms:
        grouped_tifs[chrom] = read_tifs_chrom(tifs, chrom)
    return grouped_tifs


def read_tifs_chrom(tifs, chrom, positive_strand_only=False):
    """
    For each gene in a chromosome group all tifs associated with it.

    Parameters:
        tifs (pandas): all tifs from steinmetz
        chrom (int or str): chromosome of interest
        positive_strand_only (bool): whether to filter out tifs on negative strand
    
    Returns:
        A dictionary of {gene: tifs} were tifs is a pandas with all the tifs associated with that gene
    """
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
    Extracts all reads in chrom from a bam file.

    Parameters:
        chrom (string): name of chromosome
        reads (pysam.AlignmentFile): Experiment data (bam file read in by pysam)

    Returns:
        list: List of gene objects 
    """

    filtered = reads.fetch(chrom)
    if positive_strand_only:
        return [x for x in filtered if not x.is_reverse]
    return list(reads.fetch(chrom))


def generate_5utr_isoform_starts(tifs):
    """
    Extracts all the unique 5' ends from a pandas of tifs. 
    """
    # This is a naive version for now
    return list(unique(tifs.t5))


def save_density_plot(density,
                      ox,
                      filename,
                      xlab="Nucleotides from junction",
                      ylab="Reads aligned",
                      title="Metagene"):
    """
    Saves a density plot to given filename. 
    """
    if not ox:
        ox = range(0, len(density))
    plt.plot(ox, density)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.savefig(filename)
    plt.close()


def generate_metagene(reads, tifs, chrom):
    """
    Takes all the unique TIF 5'starts from steinmetz (from genes in a given chromosome), 
    generates read densities for each gene and overlaps all the junctions into a metagene.

    Parameters:
        reads: bamfile with mRNA reads
        tifs (pandas): tif data
        chrom (string): chromosome of interest 

    Returns:
        The density used for plotted. (sum of all densities around junctions)

    Effect:
        Saves a density plot.
    """
    chrom_tifs = read_tifs_chrom(tifs, chrom, positive_strand_only=True)
    density = generate_read_density_chrom(chrom, reads)
    #TODO(astanciu): implement a funciton that generates the density of reads for a whole chromosome (no star/end params)

    splits = []

    for gene, tifs in chrom_tifs.items():
        split_positions = generate_5utr_isoform_starts(tifs)
        splits = splits + split_positions
        #list(np.random.choice(split_positions, 1, replace=False))
    return make_metagene_plot(density,
                              splits,
                              "tester_output/meta_" + chrom + ".png",
                              200,
                              200,
                              plottitle="Metagene " + chrom +
                              " steinmetz 5'UTR junction")


def generate_random_metagene(reads, chrom, sample_size):
    """
    Generates sample_size random indices inside chrom boundaries, 
    and makes a metagene around them.

    Parameters:
        reads: bamfile with mRNA reads
        chrom (string): chromosome of interest 
        sample_size (int): number of randomly generated positions

    Returns:
        The density used for plotting. (sum of all densities around junctions)

    Effect:
        Saves a density plot.
    """
    density = generate_read_density_chrom(chrom, reads)
    splits = np.random.randint(201, len(density) - 201, size=sample_size)
    # TODO(astanciu): check for duplicates
    return make_metagene_plot(density,
                              splits,
                              "tester_output/meta_" + chrom +
                              "_random_sample.png",
                              200,
                              200,
                              plottitle="Metagene " + chrom + " random sample")


def make_metagene_plot(density,
                       splits,
                       plotname,
                       upstream,
                       downstream,
                       plottitle="Metagene"):
    """ 
    Generates a metagene plot from a read density array, a list of indices and a window around each split.
    
    Sums up the read densitie in the window (-upstream, +downstream) around each split.
    


    Parameters:
        density (int list): array of read density
        splits (int list): Positions to generate metagenes around.
        plotname (string): name of plot to be saved
        upstream (int): How many nucleotides does the window include upstream of each split
        downstream (int): How many nucleotides does the window include downstream of each split 
    
    Returns:
        A meta-density array of length (downstream + upstream).

    Effect?:
        Saves a metagene plot.
    """
    print("Number of split sites is: " + str(len(splits)))
    meta = np.zeros(upstream + downstream)
    for splitidx in splits:
        meta = np.add(meta,
                      density[(splitidx - upstream):(splitidx + downstream)])

    save_density_plot(meta,
                      range(-upstream, downstream),
                      plotname,
                      title=plottitle)
    return meta


def read_files(bedfile=_kDefaultBedFile,
               bamfile=_kDefaultBamFile,
               mtiffile=_kDefaultMTIFFile,
               positive_strand_only=False):
    """
    Reads in a bedfile, a genome (bed) and steinmetz mTIF data.
    """

    # TODO(astanciu): split into multiple functions

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


def generate_read_density_chrom(chrom, allreads):
    """
    Grabs all the reads in a chrom and generates a read density array 
    (i.e. How many reads align with each position in the chrom)
    
    This function is linear in number of reads and length of chromosome.
    (doesn't depend on read length)
    """
    end = 0
    reads = grab_reads(chrom, allreads, positive_strand_only=True)
    for read in reads:
        end = max(end, read.reference_end)
    # To compute the density first compute the delta in read density at each position in the genome
    density = [0] * (end + 2)
    for read in reads:
        # The start of a read increases the delta by 1
        density[read.reference_start] = density[read.reference_start] + 1
        # the end decreases delta by 1
        density[read.reference_end + 1] = density[read.reference_end + 1] - 1
    # Summing up the deltas gives us the density.
    for i in range(1, len(density)):
        density[i] = density[i] + density[i - 1]
    return density


def generate_read_density_interval(start, end, allreads):
    """
    Grabs all the reads that overlap in [start, end] and generates a read density array 
    (i.e. How many reads align with each position in the chrom)
    
    This function is linear in number of reads and length of chromosome.
    (doesn't depend on read length)
    """
    # For more details check generate_read_density_chrom
    density = [0] * (end - start + 2)
    #TODO(astanciu): should not use chrI here but rather a parameter.
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
