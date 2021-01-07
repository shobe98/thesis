import helper
import numpy as np

np.random.seed(213)

_, genome = helper.read_files()
genome = helper.organize_genome_by_chrom(genome)
# for a gene, the start of the thick (CDS) annotated region is x.fields[6]

min_5utr_start = {}
# TODO remove this
aug_pos = {}

# For each tif, take the range [min 5'end from all tifs and annotation, annotated starting codon].
# All the values in these ranges will later used to generate random junctions
# TODO(astanciu): negative strand tifs have their ends reversed, need to compare values with annotations from genome
tifs = helper.read_all_tifs(positive_strand_only=True, return_dict=True)
for chrom, genes in tifs.items():
    min_5utr_start[chrom] = {}
    for gn, ts in genes.items():
        # for positive strand:
        min_5utr_start[chrom][gn] = min(ts['t5'])

# TODO(astanciu): add a switch for whether to include genes that weren't identified by steinmetz TIF analyisis
for chrom, genes in genome.items():
    aug_pos[chrom] = {}
    for gn, gene in genes.items():
        if gene.strand == '+':
            # Should always be the 6th field in our bed file.
            # TODO(astanciu): This might change when using a GFF.
            aug_pos[chrom][gn] = int(gene.fields[6])
            if not min_5utr_start[chrom].get(gn):
                min_5utr_start[chrom][gn] = gene.start
            else:
                if gene.start < min_5utr_start[chrom][gn]:
                    print(
                        gn +
                        " has the annotated start less than the min from TIFs "
                        + str(gene.start) + " " +
                        str(min_5utr_start[chrom][gn]))
                min_5utr_start[chrom][gn] = min(min_5utr_start[chrom][gn],
                                                gene.start)

for chrom, genes in genome.items():
    coding_length = 0
    for gn in genes.keys():
        if aug_pos[chrom].get(gn) and min_5utr_start[chrom].get(gn):
            coding_length = coding_length + (aug_pos[chrom][gn] -
                                             min_5utr_start[chrom][gn])
    print(chrom + " " + str(coding_length * 1.0 / 100000.0) +
          " total candidates for junctions")
