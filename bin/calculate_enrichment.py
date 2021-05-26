#!/bin/python
###
#   name    | mary lauren benton
#   created | 2017
#   updated | 2018.10.09
#           | 2018.10.11
#           | 2018.10.29
#           | 2019.02.01
#           | 2019.04.08
#           | 2019.06.10
#           | 2019.11.05
#           | 2021.02.24
#           | 2021.52.26
#
#   depends on:
#       BEDtools v2.23.0-20 via pybedtools
#       /dors/capra_lab/users/bentonml/data/dna/[species]/[species]_blacklist_gap.bed
#       /dors/capra_lab/data/dna/[species]/[species]-blacklist.bed
#
#       ACCRE: ml Anaconda3 GCC MySQL-client
###

import os
import sys, traceback
import argparse
import datetime
import numpy as np
from functools import partial
from multiprocessing import Pool
from pybedtools import BedTool
from pybedtools.helpers import BEDToolsError, cleanup, get_tempdir, set_tempdir


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enrichment between bed files.")

arg_parser.add_argument("region_file_1", help='bed file 1 (shuffled)')
arg_parser.add_argument("region_file_2", help='bed file 2 (not shuffled)')

arg_parser.add_argument("-i", "--iters", type=int, default=100,
                        help='number of simulation iterations; default=100')

arg_parser.add_argument("-s", "--species", type=str, default='hg19', choices=['hg19', 'hg38', 'mm10', 'dm3', 'sacCer3'],
                        help='species and assembly; default=hg19')

arg_parser.add_argument("-b", "--blacklist", type=str, default=None,
                        help='custom blacklist file; default=None')

arg_parser.add_argument("-n", "--num_threads", type=int,
                        help='number of threads; default=SLURM_CPUS_PER_TASK or 1')

arg_parser.add_argument("--print_counts_to", type=str, default=None,
                        help="print expected counts to file")

arg_parser.add_argument("--elem_wise", action='store_true', default=False,
                        help='perform element-wise overlaps; default=False')

arg_parser.add_argument("--by_hap_block", action='store_true', default=False,
                        help='perform haplotype-block overlaps; default=False')

args = arg_parser.parse_args()

# save parameters
ANNOTATION_FILENAME = args.region_file_1
TEST_FILENAME = args.region_file_2
COUNT_FILENAME = args.print_counts_to
ITERATIONS = args.iters
SPECIES = args.species
ELEMENT = args.elem_wise
HAPBLOCK= args.by_hap_block
CUSTOM_BLIST = args.blacklist

# calculate the number of threads
if args.num_threads:
    num_threads = args.num_threads
else:
    num_threads = int(os.getenv('SLURM_CPUS_PER_TASK', 1))

# if running on slurm, set tmp to runtime dir
set_tempdir(os.getenv('ACCRE_RUNTIME_DIR', get_tempdir()))


###
#   functions
###
def loadConstants(species, custom=''):
    if custom is not None:
        return custom
    return {'hg19': "/dors/capra_lab/users/bentonml/data/dna/hg19/hg19_blacklist_gap.bed",
            'hg38': "/dors/capra_lab/users/bentonml/data/dna/hg38/hg38_blacklist_gap.bed",
            'mm10': "/dors/capra_lab/users/bentonml/data/dna/mm10/mm10_blacklist_gap.bed",
            'dm3' : "/dors/capra_lab/data/dna/fly/dm3-blacklist.bed"
            }[species]


def calculateObserved(annotation, test, elementwise, hapblock):
    obs_sum = 0

    if elementwise:
        obs_sum = annotation.intersect(test, u=True).count()
    else:
        obs_intersect = annotation.intersect(test, wo=True)

        if hapblock:
            obs_sum = len(set(x[-2] for x in obs_intersect))
        else:
            for line in obs_intersect:
                obs_sum += int(line[-1])

    return obs_sum


def calculateExpected(annotation, test, elementwise, hapblock, species, custom, iters):
    BLACKLIST = loadConstants(species, custom)
    exp_sum = 0

    try:
        rand_file = annotation.shuffle(genome=species, excl=BLACKLIST, chrom=True, noOverlapping=True)

        if elementwise:
            exp_sum = rand_file.intersect(test, u=True).count()
        else:
            exp_intersect = rand_file.intersect(test, wo=True)

            if hapblock:
                exp_sum = len(set(x[-2] for x in exp_intersect))
            else:
                for line in exp_intersect:
                    exp_sum += int(line[-1])
    except BEDToolsError:
        exp_sum = -999

    return exp_sum


def calculateEmpiricalP(obs, exp_sum_list):
    mu = np.mean(exp_sum_list)
    sigma = np.std(exp_sum_list)
    dist_from_mu = [exp - mu for exp in exp_sum_list]
    p_sum = sum(1 for exp_dist in dist_from_mu if abs(exp_dist) >= abs(obs - mu))

    # add pseudocount only to avoid divide by 0 errors
    if mu == 0:
        fold_change = (obs + 1.0) / (mu + 1.0)
    else:
        fold_change = obs / mu

    p_val = (p_sum + 1.0) / (len(exp_sum_list) + 1.0)

    return "%d\t%.3f\t%.3f\t%.3f\t%.3f" % (obs, mu, sigma, fold_change, p_val)


###
#   main
###
def main(argv):
    # print header
    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())[:20]))
    print('Observed\tExpected\tStdDev\tFoldChange\tp-value')

    # run initial intersection and save
    obs_sum = calculateObserved(BedTool(ANNOTATION_FILENAME), BedTool(TEST_FILENAME), ELEMENT, HAPBLOCK)

    # create pool and run simulations in parallel
    pool = Pool(num_threads)
    partial_calcExp = partial(calculateExpected, BedTool(ANNOTATION_FILENAME), BedTool(TEST_FILENAME), ELEMENT, HAPBLOCK, SPECIES, CUSTOM_BLIST)
    exp_sum_list = pool.map(partial_calcExp, [i for i in range(ITERATIONS)])

    # wait for results to finish before calculating p-value
    pool.close()
    pool.join()

    # remove iterations that throw bedtools exceptions
    final_exp_sum_list = [x for x in exp_sum_list if x >= 0]
    exceptions = exp_sum_list.count(-999)

    # calculate empirical p value
    if exceptions / ITERATIONS <= .1:
        print(calculateEmpiricalP(obs_sum, final_exp_sum_list))
        print(f'iterations not completed: {exceptions}', file=sys.stderr)
    else:
        print(f'iterations not completed: {exceptions}\nresulted in nonzero exit status', file=sys.stderr)
        cleanup()
        sys.exit(1)

    if COUNT_FILENAME is not None:
        with open(COUNT_FILENAME, "w") as count_file:
            count_file.write('{}\n{}\n'.format(obs_sum, '\t'.join(map(str, exp_sum_list))))

    # clean up any pybedtools tmp files
    cleanup()


if __name__ == "__main__":
    main(sys.argv[1:])

