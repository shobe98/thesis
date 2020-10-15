import helper
from matplotlib import interactive
import matplotlib.pyplot as plt

interactive(True)

bamfile, bedfile, steinmetz = helper.read_files()
#a = helper.grab_reads("chrI", x, positive_strand_only=True)

#somelist = helper.generate_read_density(0, helper._kDefaultChrom1Size, bamfile)
#
#filtered = steinmetz[steinmetz.strand == "+"]
#filtered = filtered[filtered.chr == 1]
#
## first 20 genes in order of standard deviation of 5'utr length of mTIFs
#interesting = filtered.sort_values(by=["sd5"])[:20]
#
#interesting_genes = filter(lambda gene: gene.name in list(interesting.gene),
#                           bedfile)
#interesting_genes = list(interesting_genes)

#i = 0
#for gene in interesting_genes:
#    if gene.strand == "+":
#        print("Plotting gene " + gene.name)
#        plt.plot(range(1, gene.end - gene.start + 1),
#                 somelist[gene.start:gene.end])
#        plt.savefig("tester_output/gene_" + gene.name + ".png")
#        plt.close()
#
#        i = i + 1
#        if i == 20:
#            break

tifs = helper.read_all_tifs(positive_strand_only=True)
helper.generate_metagene(bamfile, bedfile, tifs, "chrI")
helper.generate_random_metagene(bamfile, bedfile, "chrI", 1231)
