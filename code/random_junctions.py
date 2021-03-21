import helper

_, genome = helper.read_files()
genome = helper.organize_genome_by_chrom(genome)

# TODO(astanciu): change to include negative strand too.
tifs = helper.read_all_tifs(positive_strand_only=True, return_dict=True)


# for a gene, the start of the thick (CDS) annotated region is x.fields[6]
def get_max_possible_genomic_regions(genome=genome, tifs=tifs):
    # FIXME assumes positive strand only.

    min_5utr_start = {}
    max_3utr_end = {}

    for chrom, genes in tifs.items():
        min_5utr_start[chrom] = {}
        max_3utr_end[chrom] = {}
        for gn, ts in genes.items():
            # for positive strand:
            min_5utr_start[chrom][gn] = min(ts['t5'])
            max_3utr_end[chrom][gn] = max(ts['t3'])

    # TODO(astanciu): add a switch for whether to include genes that weren't identified by steinmetz TIF analyisis
    for chrom, genes in genome.items():
        for gn, gene in genes.items():
            if gene.strand == '+':
                # Should always be the 6th field in our bed file.
                # TODO(astanciu): This might change when using a GFF.
                if not min_5utr_start[chrom].get(gn):
                    min_5utr_start[chrom][gn] = gene.start
                else:
                    min_5utr_start[chrom][gn] = min(min_5utr_start[chrom][gn],
                                                    gene.start)
                if not max_3utr_end[chrom].get(gn):
                    max_3utr_end[chrom][gn] = gene.end
                else:
                    max_3utr_end[chrom][gn] = max(max_3utr_end[chrom][gn],
                                                  gene.end)
    return min_5utr_start, max_3utr_end


# TODO refactor this code
def get_5utr_boundaries(genome=genome, tifs=tifs):
    min_5utr_start = {}
    # TODO remove this
    aug_pos = {}

    # For each tif, take the range [min 5'end from all tifs and annotation, annotated starting codon].
    # All the values in these ranges will later used to generate random junctions
    # TODO(astanciu): negative strand tifs have their ends reversed, need to compare values with annotations from genome
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
                    if gene.start + 1 < min_5utr_start[chrom][gn]:
                        print(
                            gn +
                            " has the annotated start less than the min from TIFs "
                            + str(gene.start) + " " +
                            str(min_5utr_start[chrom][gn]))
                    min_5utr_start[chrom][gn] = min(min_5utr_start[chrom][gn],
                                                    gene.start)

    return min_5utr_start, aug_pos


# Generates a list of all genomic postition that lie on the POSITIVE STRAND
# For each gene, consider all the positions between the minin 5'UTR end (any TIF or the annotation) and the CDS start (thick.start in the bed file)
# must take a genome organised by chromosome
def generate_candidates_for_junctions(genome=genome, tifs=tifs):
    min_5utr_start, aug_pos = get_5utr_boundaries(genome, tifs)
    candidates = {}
    for chrom, genes in genome.items():
        coding_length = 0
        cands_chrom = []
        for gn in genes.keys():
            if aug_pos[chrom].get(gn) and min_5utr_start[chrom].get(gn):
                coding_length = coding_length + (aug_pos[chrom][gn] -
                                                 min_5utr_start[chrom][gn])
                cands_chrom = cands_chrom + list(
                    range(min_5utr_start[chrom][gn], aug_pos[chrom][gn]))
        print(chrom + " " + str(coding_length * 1.0 / 100000.0) +
              " total candidates for junctions")
        candidates[chrom] = cands_chrom

    return candidates
