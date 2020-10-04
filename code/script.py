import pysam
import matplotlib.pyplot as plt
from matplotlib import interactive
from pybedtools import BedTool
import pandas as pd

interactive(True)

chr1_approx_size = 230220
# First mRNA seq file was used together with its index. Both perfectly aligned and before perfect alignment versions were used
samfile = pysam.AlignmentFile("AitkenLab/rnaseq_all.bam", "rb")
srange = list(samfile.fetch("chrI"))
somelist = [0] * (2 * chr1_approx_size)
for a in srange:
    if not a.is_reverse:
        somelist[a.reference_start] = somelist[a.reference_start] + 1
        somelist[a.reference_end + 1] = somelist[a.reference_end + 1] - 1
for i in range(1, chr1_approx_size):
    somelist[i] = somelist[i] + somelist[i - 1]

print(max(somelist))
genes = BedTool("indexes/yeast-all.bed")
stmtz = pd.read_csv("./SteinmetzGilbert/stmtz_n_mtifs.txt", sep='\t')
filtered = stmtz[stmtz.strand == "+"]
filtered = filtered[filtered.chr == 1]

# first 20 genes in order of standard deviation of 5'utr length of mTIFs
interesting = filtered.sort_values(by=["sd5"])[:20]

interesting_genes = filter(
    lambda gene: gene.name in list(interesting.gene), genes)
print(list(interesting_genes))

i = 0
for gene in interesting_genes:
    if gene.strand == "+":
        plt.plot(range(1, gene.end - gene.start + 1),
                 somelist[gene.start:gene.end])
        plt.savefig("gene_" + gene.name + ".png")
        plt.close()

        i = i + 1
        if i == 20:
            break
