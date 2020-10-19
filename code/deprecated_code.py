density = helper.generate_read_density(0, helper._kDefaultChrom1Size, bamfile)

# steinmetz is "../SteinmetzGilbert/stmtz_n_mtifs.txt"
# this code selects the 20 most interesting (highest sd of 5' length) of chr1 and plots the read distribution for each one of them
filtered = steinmetz[steinmetz.chr == 1]
# first 20 genes in order of standard deviation of 5'utr length of mTIFs
interesting = filtered.sort_values(by=["sd5"])[:20]
interesting_genes = filter(lambda gene: gene.name in list(interesting.gene),
                           bedfile.values())
interesting_genes = list(interesting_genes)

i = 0
for gene in interesting_genes:
    if gene.strand == "+":
        print("Plotting gene " + gene.name)
        plt.plot(range(1, gene.end - gene.start + 1),
                 somelist[gene.start:gene.end])
        plt.savefig("tester_output/gene_" + gene.name + ".png")
        plt.close()

        i = i + 1
        if i == 20:
            break
