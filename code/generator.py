import numpy as np
import helper
from simulate_reads import generate_isoform_aware_density, Isoform
from scipy.stats import gennorm

import pickle
import logging

#for handler in logging.root.handlers[:]:
#    logging.root.removeHandler(handler)
logging.basicConfig(filename="example.log", level=logging.WARN)

logger = logging.getLogger("logger")
logger.setLevel(logging.WARN)

shape = 1.5832492463511025
scale = 78.63922662988816

gennorm_parameters = (1.031966908880634, 6.300794747865959, 1.0705689069170423)

np.random.seed(527)

genome = helper.organize_genome_by_chrom(helper.read_bedfile())

logger.info("hello")
data = []
counter = 0
for chrom in helper.kYeastChroms:
    for gn, gene in genome[chrom].items():
        if gene.strand == "-":
            continue
        #print(gn)
        start = int(genome[chrom][gn].fields[6])
        end = int(genome[chrom][gn].fields[7])
        counts = int(np.exp(gennorm.rvs(gennorm_parameters[0],
                        loc=gennorm_parameters[1],
                        scale=gennorm_parameters[2])))
        isoforms_n = np.random.randint(1, 6)
        #print("shape ready")

        isoforms = np.random.default_rng().gamma(shape, scale, isoforms_n)

        #print("isoforms1")

        isoforms = np.unique([int(x) for x in sorted(isoforms)])

        #print("isoforms2")
        distr = [np.random.randint(1, 1000) for i in range(0, len(isoforms))]
        distr = np.divide(distr, sum(distr))

        #print("distr")

        isoforms = [Isoform(start - x, end) for x in list(isoforms)]

        #print("isoforms3")

        try:
            density = generate_isoform_aware_density(isoforms,
                                                     distr,
                                                     counts,
                                                     read_length=44,
                                                     plot=False,
                                                     crop_density=True)
        except:
            print(isoforms)
            print(distr)
            print(counts)
            raise Exception('spam', 'eggs')

        #print("density")

        data.append({
            'chrom': chrom,
            'strand': gene.strand,
            'gene': gn,
            'isoforms': isoforms,
            'distr': distr,
            'density': density
        })
        counter = counter + 1
        if counter % 50 == 0:
            print("Processed " + str(counter))
    print("Done with " + chrom)
pickle.dump(data, open("second_dataset.pickle", "wb"))
