import helper
import pickle

reads = helper.read_bamfile()
genome = helper.read_bedfile()
genome = helper.organize_genome_by_chrom(genome)


def get_bad_genes(density, chrom, strand):
    avg_readser_nt = int(sum(density) / 44.0 / len(density))
    print("average {} reads for {}".format(avg_readser_nt, chrom))
    x = [(pos, density[pos]) for pos in range(0, len(density))
         if density[pos] > avg_readser_nt * 10000]
    bad_genes = {}
    for pos, rd in x:
        for gn, gene in genome[chrom].items():
            if gene.strand == strand and gene.start <= pos and gene.end <= pos:
                if not bad_genes.get(gn):
                    bad_genes[gn] = []
                bad_genes[gn].append((pos, rd))

    return bad_genes


for chrom in helper.kYeastChroms:
    density_p, density_n = helper.generate_read_density_chrom(chrom, reads)
    bad_p = get_bad_genes(density_p, chrom, "+")
    bad_n = get_bad_genes(density_n, chrom, "-")
    print(bad_p)
    print(bad_n)
    pickle.dump({
        "+": bad_p,
        "-": bad_n
    }, open(chrom + "_bad_genes.pickle", "wb"))
