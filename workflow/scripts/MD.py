# Metagenomic diversity (MD) calculator
# Measures size and (dis)similarity of protein clusters from
# MMSeqs2 cluster and alignment outputs

# Written by Damien Finn, damien.finn@thuenen.de

import os
import sys
import json
import argparse
import numpy as np
import multiprocessing as mp
import random
from collections import defaultdict


def get_args():
    parser = argparse.ArgumentParser(
        description="Derive a Metagenomic Diversity index from clustered amino acid sequences"
    )
    parser.add_argument(
        "-c",
        "--input_clust",
        required=True,
        help="Path to mmseqs cluster results in TSV format",
    )
    parser.add_argument(
        "-m",
        "--input_aln",
        required=True,
        help="Path to mmseqs search results in .m8 format",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to write output file in JSON format",
    )
    parser.add_argument(
        "-r",
        "--random-sampling",
        action="store_true",
        help="Optional: randomly subsample observations",
    )
    parser.add_argument(
        "-N",
        "--seq-obs",
        required=False,
        type=int,
        help="""If randomly subsampling (option -r), N denotes the number of
        observations to subsample without replacement""",
    )
    parser.add_argument(
        "--random_seed",
        type=int,
        help="Random seed for random subsampling (default: None)"
    )
    parser.add_argument(
        "-T",
        "--cpus",
        default=int(mp.cpu_count()),
        type=int,
        help="Max number of CPUs (default: total available)",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    try:
        args = {
            'input_clust' : snakemake.input['clust'],
            'input_aln' : snakemake.input['aln'],
            'output' : snakemake.output[0],
            'random_sampling' : snakemake.params['random_sampling'],
            'seq_obs' : int(snakemake.params['seq_obs']),
            'random_seed' : int(snakemake.params['random_seed']),
            'cpus' : int(snakemake.threads),
        }
        sys.stderr = open(snakemake.log[0], "w") # redirect error msgs to log
    except NameError:
        args = vars(get_args())

    # Generate a dictionary of clusters
    print(" -- Making Protein dictionary -- ", file=sys.stderr)
    clustlist = []

    with open(args['input_clust']) as file:
        for l in file:
            clustlist.append(l.rstrip())

    if args['random_sampling']:
        if args['random_seed']:
            random.seed(int(args['random_seed']))
        rarvals = random.sample(clustlist, int(args['seq_obs']))
        clustlist = rarvals

    clusters = defaultdict(list)

    pool = mp.Pool(args['cpus'])


    def make_clust_dict(clustlist):
        for l in clustlist:
            col1 = l.split("\t")[0]
            col2 = l.split("\t")[1]
            clusters[col1].append(col2)


    pool.apply_async(make_clust_dict(clustlist))

    # Generate a dictionary of pairwise distances
    print(" -- Collating pair-wise dissimilarities -- ", file=sys.stderr)
    alnlist = []

    with open(args['input_aln']) as file:
        for l in file:
            alnlist.append(l.rstrip())

    pairsim = {}


    def make_pairsim_dict(alnlist):
        for l in alnlist:
            # first three columns: query, ref, id
            [col1, col2, col3] = l.split("\t")[0:3]
            # key on query-ref pairs (faster lookup than nested dict)
            pairsim[(col1, col2)] = col3


    pool.apply_async(make_pairsim_dict(alnlist))

    # Now to bring them together
    print(" -- Bringing things together -- ", file=sys.stderr)
    clustdist = defaultdict(list)


    def sorter(clusters, pairsim):
        for c in clusters:  # key name
            for a in clusters[c]:
                # look up pairwise identity values from search output
                if (c,a) in pairsim:
                    id_val = pairsim[(c,a)]  # pairwise values in q
                    clustdist[c].append(1 - float(id_val))


    pool.apply_async(sorter(clusters, pairsim))

    pool.close()
    pool.join()

    # Derive indices from the cluster dictionary

    print(" -- Deriving diversity indices -- ", file=sys.stderr)

    results = {}
    # Total Contigs
    results['total_protein_coding'] = len(clustlist)

    # Protein Richness
    results['protein_richness'] = len(clustdist)

    CDlens = [len(k) for k in clustdist.values()]
    tmp = np.sum(CDlens)

    # Shannon Diversity
    H = [(len(k) / tmp) * np.log(len(k) / tmp) for k in clustdist.values()]
    results['shannon_diversity'] = np.sum(H) * -1

    # Simpson Diversity
    S = [np.square(len(k) / tmp) for k in clustdist.values()]
    results['simpson_evenness'] = 1 - np.sum(S)

    # Metagenomic diversity
    MD = [(1 + (np.sum(k)) / len(k)) for k in clustdist.values()]
    results['protein_dissim_log10'] = np.log10(np.sum(MD))
    results['md_index'] = (1 / results['total_protein_coding']) * np.sum(MD)

    if args['output']:
        with open(args['output'], "w") as fh:
            fh.write(json.dumps(results,indent=4))
    else:
        print(json.dumps(results,indent=4), file=sys.stdout)

    print("done", file=sys.stderr)
