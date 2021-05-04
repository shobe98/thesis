import helper
import numpy as np
import pandas as pd
from random_junctions import generate_candidates_for_junctions

kJunctionCandidates = generate_candidates_for_junctions()


def select_one_random_junction(chrom, tifs):
    splits = []
    for gene, tfs in tifs[chrom].items():
        candidates = list(np.unique(sorted(tfs['t5'].tolist())))
        if len(candidates) > 1:
            splits.append((np.random.choice(candidates[1:]), "+"))
    return splits


def process_generated_data(data):
    densities = {}
    isoforms = {}
    for chrom in helper.kYeastChroms:
        end = max([
            max([iso.end for iso in row['isoforms']]) for row in data
            if row['chrom'] == chrom
        ])
        density = [0] * end
        densities[chrom] = density
        isoforms[chrom] = {}

    for row in data:
        starts = [iso.start for iso in row['isoforms']]
        ends = [iso.end for iso in row['isoforms']]
        start = min(starts)
        end = max(ends)
        isoforms[row['chrom']][row['gene']] = pd.DataFrame(data={
            't5': starts,
            't3': ends
        })
        densities[row['chrom']][start:(end + 2)] = np.add(
            densities[row['chrom']][start:(end + 2)], row['density'])
    return densities, isoforms
