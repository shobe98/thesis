import helper
import logging

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
else:
    logging.basicConfig(level=logging.WARN)

_, genome = helper.read_files()
genome = helper.organize_genome_by_chrom(genome)

# generate some statistics about the bam files !!!!!
# 1. how many intergenic reads do we have?
# 2. is there an increase in reads in CDS?
# 3. average read length


# TODO(astanciu): move this function to a different file
def generate_read_density_chrom(chrom, allreads):
    """
    Grabs all the reads in a chrom and generates a read density array 
    (i.e. How many reads align with each position in the chrom)
    
    This function is linear in number of reads and length of chromosome.
    (doesn't depend on read length)
    """
    logging.info("generate_read_density_chrom chrom=" + chrom)
    end = 0
    reads = grab_reads(chrom, allreads)
    for read in reads:
        end = max(end, read.reference_end)
    if DEBUG_MODE:
        logging.info("len(reads)=" + str(len(reads)))
        logging.info("chrom_size=" + str(end))

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
    logging.info("Done generating read density...")
    logging.info("Sum of density_n: " + str(sum(density_n)))
    logging.info("Sum of density_p: " + str(sum(density_p)))
    return (density_p, density_n)
