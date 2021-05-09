import helper
from random_junctions import get_max_possible_genomic_regions, get_5utr_boundaries
import matplotlib.pyplot as plt
import numpy as np

tifs = helper.read_all_tifs(positive_strand_only=False, return_dict=True)
genome = helper.organize_genome_by_chrom(helper.read_bedfile())
bamfile = helper.read_bamfile()

#gene_start, gene_end = get_5utr_boundaries()
gene_start, gene_end = get_max_possible_genomic_regions()

for chrom in helper.kYeastChroms:
    density, _ = helper.generate_read_density_chrom(chrom, bamfile)
    del _
    for gn, gene in genome[chrom].items():
        if gene.strand == "+" and gn in tifs[chrom] and gn in gene_start[
                chrom] and gn in gene_end[chrom]:
            print(gn)
            plt.plot(range(gene_start[chrom][gn], gene_end[chrom][gn]),
                     density[gene_start[chrom][gn]:gene_end[chrom][gn]])
            junctions = np.unique(sorted(tifs[chrom][gn]['t5'].tolist()))
            for j in junctions:
                plt.axvline(j, color='k', alpha=0.5)

            plt.xlabel("Genomic position (nt)")
            plt.ylabel("Read density")

            plt.title("Read density across 5'UTR of gene " + gn)

            plt.savefig("bulk/gene_" + gn + ".png")
            plt.close()
