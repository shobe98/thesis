import pickle 
import helper

reads = helper.read_bamfile()

densities = {}
for chrom in helper.kYeastChroms:
        densities[chrom] = helper.generate_read_density_chrom(chrom, reads)

f = open("densities.pickle", "wb")
pickle.dump(densities, f)
f.close()
