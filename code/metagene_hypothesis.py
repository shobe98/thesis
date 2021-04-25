import logging
import argparse
import csv
import helper

parser = argparse.ArgumentParser()

me = parser.add_subparsers(dest='input_type')
sim = me.add_parser("simulated",
                    help="Parses a pickle file of simulated reads")
sim.add_argument("input_pickle", help="Input path to pickle file.")

rna = me.add_parser("rnaseq", help="Parses a bam file of RNASeq reads")
rna.add_argument("input_bam", help="Path to bam file")
rna.add_argument("input_mtif", help="Path to annotated tifs.")

parser.add_argument("output", help="File to output dataset.")

parser.add_argument(
    "-n",
    type=int,
    default=100,
    help="Number of positions in the metagene to be considered.")
parser.add_argument("--iterations",
                    type=int,
                    default=10000,
                    help="Number of iterations for random simulations.")

args = parser.parse_args()

output_file = open(args.output, mode='w')
output_writer = csv.writer(outpur_file,
                           delimiter=',',
                           quotechar='"',
                           quoting=csv.QUOTE_MINIMAL)

if args.input_type == "rnaseq":
    tifs = helper.read_all_tifs(args.input_mtif,
                                positive_strand_only=True,
                                return_dict=False)
    reads = helper.read_bamfile(args.input_bam)

    densities = []
    for chrom in helper.kYeastChroms:
        density, _ = helper.generate_read_density_chrom(chrom, reads)
        densities.append(density)

    del _
else:
    pass

for it in range(0, args.N):
    for chrom in helper.kYeastChroms:
        pass

output_file.close()
