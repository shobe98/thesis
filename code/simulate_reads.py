import numpy as np
import helper
import random_junctions

import seaborn as sns
import matplotlib.pyplot as plt

np.random.seed(527)


class FakeRead:
    def __init__(self, start, end, strand="+"):
        self.reference_start = start
        self.reference_end = end
        self.is_reverse = (strand == "-")


def random_reads_in_interval(lo, hi, sample_size, read_length, weights=None):
    if weights is None:
        # TODO change this to choice
        reads = np.random.randint(lo, hi - read_length, sample_size)
    else:
        p = np.divide(weights[lo:hi], sum(weights[lo:hi]))
        reads = np.random.choice(range(lo, hi), size=sample_size, p=p)
    reads = [FakeRead(x, x + read_length) for x in reads]
    reads.append(FakeRead(
        hi, hi))  # quick hack so the read densities have the same size
    # FIXME
    return reads


def generate_random_reads(genome=helper.organize_genome_by_chrom(
    helper.read_bedfile()),
                          allreads=helper.read_bamfile()):
    """ for each gene generates uniformly distributed reads """
    density_p_original, _ = helper.generate_read_density_chrom(
        "chrI", allreads)
    for gname, gene in genome["chrI"].items():
        if gene.strand == "+":
            sample_size = int(
                sum(density_p_original[gene.start:gene.end]) / 44)
            fake, _ = helper.__generate_read_density(
                random_reads_in_interval(gene.start, gene.end, sample_size,
                                         44))

            fake10, _ = helper.__generate_read_density(
                random_reads_in_interval(gene.start, gene.end,
                                         sample_size * 10, 44))
            plt.plot(range(gene.start, gene.end),
                     density_p_original[gene.start:gene.end], "b")
            plt.plot(range(gene.start, gene.end), fake[gene.start:-2], "r")
            plt.plot(range(gene.start, gene.end), fake10[gene.start:-2], "g")
            plt.title(
                gname +
                " Real data (b) vs randomly generated reads, and 10x randomly generated reads"
            )
            plt.savefig("random_generation/individual_genes/" + gname + ".png")
            plt.close()


def generate_uniform_reads_for_chrom():
    """ for each chromosome, generates uniformly distributed reads. Uses the original reads and the original reads with filtered outliers """
    allreads = helper.read_bamfile()
    for chrom in helper.kYeastChroms:
        print("Generating density from bam")

        density_p, _ = helper.generate_read_density_chrom(chrom, allreads)
        __generate_uniform_reads_for_chrom(
            chrom, density_p, "random_generation/end_first_no_filter")

        # Attempting to remove outliers
        density_p = remove_outliers(density_p)
        __generate_uniform_reads_for_chrom(
            chrom, density_p, "random_generation/end_first_filter")

        print("stopping after chrI... please disable this")
        break


def generate_weighted_reads():
    allreads = helper.read_bamfile()
    for chrom in helper.kYeastChroms:
        print("Generating density from bam")

        density_p, _ = helper.generate_read_density_chrom(chrom, allreads)
        weights = generate_weights_smart(density_p, chrom)


def fix_lengths(a, b):
    if len(a) < len(b):
        a = a + [0] * (len(b) - len(a))
    if len(b) < len(a):
        b = b + [0] * (len(a) - len(b))
    return a, b


def __generate_uniform_reads_for_chrom(chrom, density_p_original,
                                       plot_location):
    sample_size = int(sum(density_p_original) / 44)
    print("Generating %d reads..." % sample_size)
    print("Generating uniform data...")
    uniform, _ = helper.__generate_read_density(
        random_reads_in_interval(0, len(density_p_original), sample_size, 44))
    uniform, density_p_original = fix_lengths(uniform, density_p_original)

    print("Generating neively weighted data...")
    weighted_naive, _ = helper.__generate_read_density(
        random_reads_in_interval(0, len(density_p_original), sample_size, 44,
                                 generate_weights_naive(density_p_original)))
    weighted_naive, density_p_original = fix_lengths(weighted_naive,
                                                     density_p_original)

    print("Generating weighted data (using the avg method)...")
    weighted_avg, _ = helper.__generate_read_density(
        random_reads_in_interval(
            0, len(density_p_original), sample_size, 44,
            generate_weights_averaged(density_p_original)))
    weighted_avg, density_p_original = fix_lengths(weighted_avg,
                                                   density_p_original)

    for i in range(0, len(uniform), 2000):
        print("Plotting from " + str(i) + "...")
        end = min(i + 2000, len(uniform))
        plt.plot(range(i, end),
                 density_p_original[i:end],
                 "b",
                 label="original density")
        plt.plot(range(i, end),
                 uniform[i:end],
                 "r",
                 label="uniformly generated density")
        plt.plot(range(i, end),
                 weighted_naive[i:end],
                 "g",
                 label="naively weighted density")
        plt.plot(range(i, end),
                 weighted_avg[i:end],
                 "y",
                 label="averaged out weights - density")
        plt.title("Genome wide \n Real data (b) vs randomly generated reads")
        plt.savefig(plot_location + "/{0}_density_line_".format(chrom) +
                    str(i) + ".png")
        plt.close()

    # TODO move to separate function
    #print("Making dataset for histogram...")
    #print("og: %d, uniform: %d" % (len(density_p_original), len(uniform)))
    ##data = list(np.subtract(density_p_original, uniform))
    #data = list(np.log(density_p_original))
    #print(len(data))
    #if len(data) % 10000 != 0:
    #    data = data + [0] * (10000 - len(data) % 10000)
    #print(len(data))
    #data = [data[i:(i + 10000)] for i in range(0, len(data), 10000)]
    #print(len(data))

    #sns.heatmap(data).figure.savefig(
    #    "random_generation/whole_genome/{0}_heatmap_of_our_data_log_scale.png".
    #    format(chrom))
    #plt.close()
    print("Done")


def generate_weights_naive(density):
    print("Generating naive weights...")
    # each position will be selected as a START for
    # an interval proportional with the number of reads aligned to that position
    weights = np.divide(density, sum(density))
    return weights


def generate_weights_averaged(density):
    print("Generating avg weights...")
    # each position will be selected as a START for an interval
    # proportional to the average number of reads starting at a position
    # The average number of reads is the number of reads aligning divided by the avg read length
    weights = np.divide(density, sum(density) / 44)
    return weights


def remove_outliers(density):
    outlier_thr = np.quantile(density, 0.99)
    # 1 if good, 0 if bad
    outliers = [
        int(density[pos] < outlier_thr) for pos in range(0, len(density))
    ]
    # keeps the original values for good positions and completely erases the bad ones.
    density = list(np.multiply(density, outliers))
    return density


def generate_weights_smart(density,
                           chrom,
                           genome=helper.organize_genome_by_chrom(
                               helper.read_bedfile())):
    weights = remove_outliers(density)
    pos_5p, pos_3p = random_junctions.get_max_possible_genomic_regions(genome)

    genic = [0] * len(density)
    for gn, gene in genome[chrom].items():
        if gene.strand == "+":
            if pos_5p.get(gn):
                r = range(pos_5p[gn], pos_3p[gn] + 1)
            else:
                r = range(gene.start, gene.end + 1)
            for i in r:
                genic[i] = 1

    genic_reads = list(np.multiply(density, genic))
    print("Genic {}/{}%, intergenic {}, total {}".format(
        sum(genic_reads) / 44.0,
        sum(genic_reads) / sum(density),
        (sum(density) - sum(genic_reads)) / 44, sum(density)))
    return None


#generate_random_reads()
generate_uniform_reads_for_chrom()
#generate_weighted_reads()
