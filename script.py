import pysam
import matplotlib.pyplot as plt
from matplotlib import interactive
from pybedtools import BedTool

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

i = 0
for agene in genes:
    if agene.strand == "+":
        plt.plot(range(1, agene.end - agene.start + 1),
                 somelist[agene.start:agene.end])
        plt.savefig("gene" + str(i) + ".png")
        plt.close()

        i = i + 1
        if i == 20:
            break
