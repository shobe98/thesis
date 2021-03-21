import numpy as np
import helper

import seaborn as sns
import matplotlib.pyplot as plt

np.random.seed(527)


class FakeRead:
    def __init__(self, start, end, strand="+"):
        self.reference_start = start
        self.reference_end = end
        self.is_reverse = (strand == "-")


def random_reads_in_interval(lo, hi, sample_size, read_length):
    reads = np.random.randint(lo, hi - read_length, sample_size)
    reads = [FakeRead(x, x + read_length) for x in reads]
    reads.append(FakeRead(
        hi, hi))  # quick hack so the read densities have the same size
    # FIXME
    return reads


def generate_random_reads(genome=helper.organize_genome_by_chrom(
    helper.read_bedfile()),
                          allreads=helper.read_bamfile()):
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
    for chrom in helper.kYeastChroms:
        __generate_uniform_reads_for_chrom(chrom)


def __generate_uniform_reads_for_chrom(chrom, allreads=helper.read_bamfile()):
    print("Generating density from bam")
    density_p_original, _ = helper.generate_read_density_chrom(chrom, allreads)
    sample_size = int(sum(density_p_original) / 44)
    print("Generating %d reads..." % sample_size)
    print("Generating fake data...")
    fake, _ = helper.__generate_read_density(
        random_reads_in_interval(0, len(density_p_original), sample_size, 44))
    if len(fake) < len(density_p_original):
        fake = fake + [0] * (len(density_p_original) - len(fake))
    if len(density_p_original) < len(fake):
        density_p_original = density_p_original + [0] * (
            len(fake) - len(density_p_original))

    for i in range(0, len(fake), 10000):
        print("Plotting from " + str(i) + "...")
        end = min(i + 10000, len(fake))
        plt.plot(range(i, end), density_p_original[i:end], "b")
        plt.plot(range(i, end), fake[i:end], "r")
        plt.title("Genome wide \n Real data (b) vs randomly generated reads")
        plt.savefig(
            "random_generation/whole_genome/{0}_density_line_".format(chrom) +
            str(i) + ".png")
        plt.close()

    print("Making dataset for histogram...")
    print("og: %d, fake: %d" % (len(density_p_original), len(fake)))
    #data = list(np.subtract(density_p_original, fake))
    data = list(np.log(density_p_original))
    print(len(data))
    if len(data) % 10000 != 0:
        data = data + [0] * (10000 - len(data) % 10000)
    print(len(data))
    data = [data[i:(i + 10000)] for i in range(0, len(data), 10000)]
    print(len(data))

    sns.heatmap(data).figure.savefig(
        "random_generation/whole_genome/{0}_heatmap_of_our_data_log_scale.png".
        format(chrom))

    print("Done")


#generate_random_reads()
generate_uniform_reads_for_chrom()
