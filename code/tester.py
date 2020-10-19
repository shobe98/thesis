import helper
from matplotlib import interactive
import matplotlib.pyplot as plt
import numpy as np

interactive(True)

_kDefaultChrom1Size = 229424

# Test for read_files (no assertions; if it works, it works :)) )
bamfile, genes_dict, steinmetz = helper.read_files()

# Test for grab_reads
chr1_reads = helper.grab_reads("chrI", bamfile, positive_strand_only=True)

# Test for generate_read_density
density = helper.generate_read_density_interval(0, _kDefaultChrom1Size,
                                                bamfile)
density2 = helper.generate_read_density_chrom("chrI", bamfile)
assert np.equal(density, density2).all()

# Test for read all tifs
tifs = helper.read_all_tifs(positive_strand_only=True)

# Test for generate metagene:
# Overlaps all junctions found in the tif file and creates a metagene plot for -200 +200 nt of each junction
# tests generate_5utr_isoform_starts and make_metagene_plot
helper.generate_metagene(bamfile, tifs, "chrI")
helper.generate_random_metagene(bamfile, "chrI", 1231)
