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
_kMetageneRangeNt = 50

DEBUG_MODE = True

kYeastChroms = [
    "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
    "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"
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
        print("Done reading TIF file... returning a pandas")
        return tifs

    grouped_tifs = {}
    for chrom in kYeastChroms:
        grouped_tifs[chrom] = read_tifs_chrom(tifs, chrom)
    print("Done reading TIF file... returning a dictionary")
    return grouped_tifs


def organize_genome_by_chrom(genome):
    grouped = {}
    for c in kYeastChroms:
        grouped[c] = {}

    for gene_id, gene in genome.items():
        grouped[gene.chrom][gene_id] = gene

    return grouped


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
    print("Processed TIFs in chromosome " + str(chrom) + "...")
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
    chr_reads = list(reads.fetch(chrom))
    if DEBUG_MODE:
        print("Processed mRNAs from chromosome " + chrom +
              ". Total reads in chromosome: " + str(len(chr_reads)) +
              ". positive_strand_only=" + str(positive_strand_only))
    return list(filtered)  # TODO(astanciu): inefficient?


def generate_5utr_isoform_starts(tifs):
    """
    Extracts all the unique 5' ends from a pandas of tifs. 
    """

    ret_p = []
    ret_n = []
    for i, x in tifs.iterrows():
        if x['strand'] == "+":
            ret_p.append(x['t5'])
        else:
            ret_n.append(x['t3'])

    # TODO positive_strand_only
    ret = [(x, "+") for x in sorted(list(unique(ret_p)))
           ]  #+ [(x, "-") for x in sorted(list(unique(ret_n)))]
    # This is a naive version for now
    return ret
    #if len(ret) <= 2:
    #    return []  # never return the first junction site
    #else:
    #    return [ret[int(len(ret) / 2)]]


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

    print("Saved plot " + filename)


def generate_metagene(reads, tifs, chrom, density=None):
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
    print("Generating metagene for chromosome " + chrom)
    chrom_tifs = read_tifs_chrom(tifs, chrom)

    if not density:
        density = generate_read_density_chrom(chrom, reads)
    print("Done generating density for metagene plot for chromosome " + chrom)
    print(type(density))
    splits = []

    for gene, tifs in chrom_tifs.items():
        split_positions = generate_5utr_isoform_starts(tifs)
        splits = splits + split_positions
        if not split_positions:
            print("WTFFF...")
            print(gene)
            print(tifs)
        #list(np.random.choice(split_positions, 1, replace=False))
    print("type of splits[0] in make_metagene " + str(type(splits[0])))
    print(splits[0])
    print("Making the metagene plot...")
    return make_metagene_plot(
        density,
        splits,
        "tester_output_w_junctions/meta_pos_" + chrom + ".png",
        _kMetageneRangeNt,
        _kMetageneRangeNt,
        plottitle="Metagene " + chrom + " steinmetz 5'UTR junction")


def inside_annotated_regions(x, strand, genome_chrom):
    return x > 200


def generate_random_split_junctions(size, candidates, chrom_size,
                                    genome_chrom):
    print("generating random junctions...")
    splits_p = np.random.choice(candidates, size, replace=False)
    splits_p = [(x, "+") for x in splits_p
                if inside_annotated_regions(x, "+", genome_chrom)]
    splits_n = []
    splits_n = list(unique(np.random.randint(1, chrom_size - 201, size)))
    splits_n = [(x, "-") for x in splits_n
                if inside_annotated_regions(x, "+", genome_chrom)]

    # TODO positive_strand_only
    return splits_p
    return splits_n + splits_p


def generate_random_metagene(reads,
                             chrom,
                             sample_size,
                             candidates,
                             density=None,
                             genome=None):
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
    print("Generating metagene for chromosome " + chrom +
          " with randomly sampled junction sites")
    if not density:
        density = generate_read_density_chrom(chrom, reads)
    if genome:
        splits = generate_random_split_junctions(size=sample_size,
                                                 candidates=candidates,
                                                 chrom_size=len(density[0]),
                                                 genome_chrom=genome[chrom])
    else:
        splits_p = list(
            np.random.randint(201, chrom_size - 201, size=sample_size))
        print(splits_p)
        splits_p = [(x, "+") for x in splits_p]
        print(splits_p)
        splits_n = list(
            np.random.randint(201, chrom_size - 201, size=sample_size))
        splits_n = [(x, "-") for x in splits_n]

        splits = splits_p + splits_n

    print(
        "Done generating density for random-sample-metagene plot for chromosome "
        + chrom)
    print(splits)
    # TODO(astanciu): check for duplicates
    print("Generating random sample metagene...")
    return make_metagene_plot(density,
                              splits,
                              "tester_output_w_junctions/meta_pos_" + chrom +
                              "_random_sample.png",
                              _kMetageneRangeNt,
                              _kMetageneRangeNt,
                              plottitle="Metagene " + chrom + " random sample")


def make_metagene_plot(density,
                       splits,
                       plotname,
                       upstream,
                       downstream,
                       plottitle="Metagene"):
    """ 
    Generates a metagene plot from a read density array, a list of indices and a window around each split.
    
    Sums up the read density in the window (-upstream, +downstream) around each split.
    


    Parameters:
        density (int list): array of read density
        splits (int * str list): positions and their corresponding strand, used to generate metagenes around.
        plotname (string): name of plot to be saved
        upstream (int): How many nucleotides does the window include upstream of each split
        downstream (int): How many nucleotides does the window include downstream of each split 
    
    Returns:
        A meta-density array of length (downstream + upstream).

    Effect?:
        Saves a metagene plot.
    """
    print("Generating metagene plot... Number of overlapping sites is: " +
          str(len(splits)))
    meta = np.zeros(upstream + downstream)
    print(type(density))
    density_p, density_n = density

    print("type of splits[0] is " + str(type(splits[0])))
    print(splits[0])
    for splitidx, strand in splits:
        if strand == "+":
            meta = np.add(
                meta, density_p[(splitidx - upstream):(splitidx + downstream)])
        else:
            meta = np.add(
                meta, density_n[(splitidx - upstream):(splitidx + downstream)])

    print("Saving metagene plot...")
    save_density_plot(meta,
                      range(-upstream, downstream),
                      plotname,
                      title=plottitle)
    return meta


def read_bamfile(bamfile=_kDefaultBamFile):
    print("Reading bamfile...")
    # This is just an iterator in the file
    bam = pysam.AlignmentFile(bamfile)
    print("Done reading bamfile " + bamfile)
    return bam


def read_bedfile(bedfile=_kDefaultBedFile):
    print("Reading bedfile...")
    genes = BedTool(bedfile)  # Full genome loaded in memory
    genes_dict = {}
    for g in genes:
        genes_dict[g.name] = g
    print("Done reading bedfile " + bedfile)
    return genes_dict


def read_files(bedfile=_kDefaultBedFile, bamfile=_kDefaultBamFile):
    """
    Reads in a bedfile, a genome (bed) and steinmetz mTIF data.
    """
    return read_bamfile(bamfile), read_bedfile(bedfile)


def read_mtifs(mtiffile=_kDefaultMTIFFile, positive_strand_only=False):
    """
    Reads in steinmetz mTIF data.
    """
    mtifs = pd.read_csv(mtiffile,
                        sep='\t')  # Full pandas dataset loaded in memory

    if positive_strand_only:
        mtifs = mtifs[mtifs.strand == "+"]
    return mtifs


def generate_read_density_chrom(chrom, allreads):
    """
    Grabs all the reads in a chrom and generates a read density array 
    (i.e. How many reads align with each position in the chrom)
    
    This function is linear in number of reads and length of chromosome.
    (doesn't depend on read length)
    """
    print("generate_read_density_chrom chrom=" + chrom)
    reads = grab_reads(chrom, allreads)
    return __generate_read_density(reads)


def __generate_read_density(reads):
    end = 0
    for read in reads:
        end = max(end, read.reference_end)
    if DEBUG_MODE:
        print("len(reads)=" + str(len(reads)))
        print("chrom_size=" + str(end))

    # To compute the density first compute the delta in read density at each position in the genome
    density_p = [0] * (end + 2)
    density_n = [0] * (end + 2)
    for read in reads:
        if read.is_reverse:
            # Reads on negative strand have start < end. I'll leave the if-else here for more fine grained control
            # The start of a read increases the delta by 1
            density_n[
                read.reference_start] = density_n[read.reference_start] + 1
            # the end decreases delta by 1
            density_n[read.reference_end +
                      1] = density_n[read.reference_end + 1] - 1
        else:
            # The start of a read increases the delta by 1
            density_p[
                read.reference_start] = density_p[read.reference_start] + 1
            # the end decreases delta by 1
            density_p[read.reference_end +
                      1] = density_p[read.reference_end + 1] - 1
    # Summing up the deltas gives us the density.
    for i in range(1, len(density_p)):
        density_p[i] = density_p[i] + density_p[i - 1]
        density_n[i] = density_n[i] + density_n[i - 1]
    print("Done generating read density...")
    print("Sum of density_n: " + str(sum(density_n)))
    print("Sum of density_p: " + str(sum(density_p)))
    return (density_p, density_n)


def generate_read_density_interval(start, end, allreads):
    """
    DEPRECATED!!!!
    Grabs all the reads that overlap in [start, end] and generates a read density array 
    (i.e. How many reads align with each position in the chrom)
    
    This function is linear in number of reads and length of chromosome.
    (doesn't depend on read length)
    """
    # For more details check generate_read_density_chrom
    print("This function is not used outside of its own test")
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
