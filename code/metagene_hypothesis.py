import logging
import argparse
import csv
import helper
import metagene_helper
import pickle
import numpy as np
from metagene_helper import select_one_random_junction, kJunctionCandidates

logger = logging.getLogger("root")
logger.setLevel(logging.INFO)

np.random.seed(100612)

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

print(args)

output_file = open(args.output, mode='w')
output_writer = csv.writer(output_file,
                           delimiter=',',
                           quotechar='"',
                           quoting=csv.QUOTE_MINIMAL)

output_writer.writerow(["iteration"] +
                       ["j_delta_" + str(i) for i in range(1, args.n)] +
                       ["f_delta_" + str(i) for i in range(1, args.n)])

if args.input_type == "rnaseq":
    tifs = helper.read_all_tifs(args.input_mtif,
                                positive_strand_only=True,
                                return_dict=True)
    reads = helper.read_bamfile(args.input_bam)

    densities = {}
    for chrom in helper.kYeastChroms:
        density, _ = helper.generate_read_density_chrom(chrom, reads)
        densities[chrom] = density
        break
    del _
else:
    generated = pickle.load(open(args.input_pickle, "rb"))
    densities, tifs = metagene_helper.process_generated_data(generated)

for it in range(0, args.iterations):
    big_meta = [0] * (2 * args.n)
    fake_meta = [0] * (2 * args.n)

    for chrom in helper.kYeastChroms:
        splits = select_one_random_junction(chrom, tifs)
        meta = helper.make_metagene_plot((densities[chrom], []),
                                         splits,
                                         "",
                                         args.n,
                                         args.n,
                                         plot=False)
        big_meta = np.add(big_meta, meta)

        fake_junctions = [(x, "+") for x in np.random.choice(
            metagene_helper.kJunctionCandidates[chrom],
            size=len(splits),
            replace=False)]
        meta = helper.make_metagene_plot((densities[chrom], []),
                                         fake_junctions,
                                         "",
                                         args.n,
                                         args.n,
                                         plot=False)
        fake_meta = np.add(fake_meta, meta)

        break
    data = [
        big_meta[args.n + i] / big_meta[args.n - i] for i in range(1, args.n)
    ] + [
        fake_meta[args.n + i] / fake_meta[args.n - i]
        for i in range(1, args.n)
    ]

    output_writer.writerow([it] + data)

output_file.close()
