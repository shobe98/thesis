import helper
import random_junctions
from pandas import Series

genome = helper.read_bedfile()
genome = helper.organize_genome_by_chrom(genome)
allreads = helper.read_bamfile()
tifs = random_junctions.tifs
# generate some statistics about the bam files !!!!!
# 1. how many intergenic reads do we have?
# 2. is there an increase in reads in CDS?
# 3. average read length


def analyze_chromosome(chrom, density_p, annotated_regions, window_size=100):
    # FIXME only works for positive strand
    genic = [0] * len(density_p)
    start, end = annotated_regions
    # mark the beginning of each read
    for gn, st in start[chrom].items():
        genic[st] = genic[st] + 1

    rolling_avg_reads = Series(density_p).rolling(window_size).mean().tolist()
    for gn, st in start[chrom].items():
        en = end[chrom][gn]
        print(gn)
        print(rolling_avg_reads[(st - 400):(st - 200)])
        print(rolling_avg_reads[st:(st + 200)])

        print(sum(density_p[st:en]))
        print("\n")


def bam_stats(allreads, genome, tifs):
    print("computing genomic boundaries")
    min_5utr_start, aug_pos = random_junctions.get_5utr_boundaries(
        genome, tifs)
    annotated_regions = random_junctions.get_max_possible_genomic_regions(
        genome, tifs)

    result = {}
    for chrom in helper.kYeastChroms:
        print("processing chrom " + chrom)
        result[chrom] = {}
        end = 0

        reads = helper.grab_reads(chrom, allreads, positive_strand_only=True)
        total_size = 0.0

        print("computing avg read size")
        for read in reads:
            end = max(end, read.reference_end)
            total_size = total_size + read.reference_length

        average_size = total_size / len(reads)

        result[chrom]["avg_size"] = average_size
        result[chrom]["total_size"] = total_size
        result[chrom]["read_count"] = len(reads)

        print("computing density of reads")
        # To compute the density first compute the delta in read density at each position in the genome
        density_p = [0] * (end + 2)
        density_n = [0] * (end + 2)
        for read in reads:
            if read.is_reverse:
                # Reads on negative strand have start < end. I'll leave the if-else here for more fine grained control
                # The start of a read increases the delta by 1
                density_n[
                    read.reference_start] = density_n[read.reference_start] + 1
            else:
                # The start of a read increases the delta by 1
                density_p[
                    read.reference_start] = density_p[read.reference_start] + 1

        print("analyze_chromosome")
        # need a "density" array just marks with 1 the beginning of any read.
        analyze_chromosome(chrom, density_p, annotated_regions)
        #this is necessary to answer queries of the form "how many reads start in given interval"
        for i in range(1, len(density_n)):
            density_p[i] = density_p[i - 1] + density_p[i]
            density_n[i] = density_n[i - 1] + density_n[i]

        print("computing metrics")
        in_annotated_regions = 0
        in_cds = 0
        in_5utr = 0
        right_before_utr = 0
        for gn, gene in genome[chrom].items():
            if gene.strand == "+":
                in_annotated_regions = in_annotated_regions + density_p[
                    gene.end] - density_p[gene.start - 1]
                in_cds = in_cds  #TODO
                if min_5utr_start[chrom].get(gn):
                    in_5utr = in_5utr + density_p[aug_pos[chrom][
                        gn]] - density_p[min_5utr_start[chrom][gn] - 1]
                    right_before_utr = right_before_utr + density_p[
                        min_5utr_start[chrom][gn]] - density_p[
                            min_5utr_start[chrom][gn] - 201]
                else:
                    print("gene " + gn + " doesn't have any info on 5'utr???")

            else:
                pass

        result[chrom]["in_annotated_regions"] = in_annotated_regions
        result[chrom]["in_cds"] = in_cds
        result[chrom]["in_5utr"] = in_5utr
        result[chrom]["right_before_utr"] = right_before_utr
        print(sum(density_n))
        break
    return result


bam_stats(allreads, genome, tifs)
