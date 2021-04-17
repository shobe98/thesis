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


class Isoform:
    # probably more data
    def __init__(self, start, end, strand="+"):
        self.start = start
        self.end = end
        self.strand = strand


# JUST PROTOTYPES
def generate_isoform_aware_density(isoforms,
                                   distribution,
                                   total_reads,
                                   read_length=44):
    if len(isoforms) != len(distribution):
        print("Error distribution and isoforms should have the same size")
        return
    isoform_reads = np.multiply(distribution, total_reads)
    isoform_reads = [int(x) for x in isoform_reads]
    min_p = min([iso.start for iso in isoforms])
    max_p = max([iso.end for iso in isoforms])
    density = [0] * (max_p - min_p + 1)
    for i in range(0, len(isoforms)):
        reads = random_reads_in_interval(isoforms[i].start, isoforms[i].end,
                                         isoform_reads[i], read_length)
        dens, _ = helper.__generate_read_density(reads)
        print(len(dens))
        d_start = isoforms[i].start - min_p
        d_end = isoforms[i].end - min_p + 1
        print(d_start)
        print(d_end)
        print(isoforms[i].start)
        print(len(dens[isoforms[i].start:-1]))
        density[d_start:d_end] = list(
            np.add(density[d_start:d_end], dens[isoforms[i].start:-1]))

    plt.plot(range(min_p, max_p + 1), density)
    plt.title("Isoform distribution\n{}".format(str(distribution)))
    plt.savefig("tinkering/isoforms/{}_{}.png".format(len(isoforms),
                                                      total_reads))
    plt.close()


def random_reads_by_fragmentation(lo, hi, sample_size, read_length):
    gene = list(range(lo, hi))
    total = []
    for i in range(0, sample_size):
        fragments = np.random.choice(gene, int(len(gene) / read_length))
        fragments = [0] + sorted(fragments) + [len(gene)]
        lengths = [
            fragments[j] - fragments[j - 1] for j in range(1, len(fragments))
        ]

        total = total + lengths

        plt.hist(total, bins=200)
        plt.savefig("tinkering/fragments_2000/" + str(i) + ".png")
        plt.close()


# ACTUAL CODE
def random_reads_in_interval(lo, hi, sample_size, read_length, weights=None):
    if weights is None:
        # TODO change this to choice
        reads = np.random.randint(lo + read_length, hi, sample_size)
    else:
        p = np.divide(weights[(lo + read_length):hi],
                      sum(weights[(lo + read_length):hi]))
        reads = np.random.choice(range(lo + read_length, hi),
                                 size=sample_size,
                                 p=p)
    reads = [FakeRead(x - read_length, x) for x in reads]
    reads.append(FakeRead(lo, lo))
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


if __name__ == "__main__":
    #generate_random_reads()
    #generate_uniform_reads_for_chrom()
    #generate_weighted_reads()
    #random_reads_by_fragmentation(0, 2000, 200, 44)

    isoforms2 = [Isoform(1000, 2000), Isoform(1200, 2000)]
    prob2 = [0.4, 0.6]
    isoforms3 = [Isoform(1000, 2000), Isoform(1200, 2000), Isoform(1250, 2000)]
    prob3 = [0.3, 0.3, 0.4]
    isoforms4 = [
        Isoform(1000, 2000),
        Isoform(1200, 2000),
        Isoform(1250, 2000),
        Isoform(1400, 2000)
    ]
    prob4 = [0.1, 0.1, 0.4, 0.4]

    for total_reads in range(100, 10000, 100):
        generate_isoform_aware_density(isoforms2, prob2, total_reads)
        generate_isoform_aware_density(isoforms3, prob3, total_reads)
        generate_isoform_aware_density(isoforms4, prob4, total_reads)
