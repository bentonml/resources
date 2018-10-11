#!/bin/python
###
#   name    | mary lauren benton
#   created | 2017
#   updated | 2018.10.09
#           | 2018.10.11
#
#   depends on:
#       BEDtools v2.23.0-20
#       /dors/capra_lab/data/dna/[species]_trim.chrom.sizes
#       /dors/capra_lab/users/bentonml/data/[species]_blacklist_gap.bed
###

import os
import sys, traceback
import argparse
import datetime
import numpy as np
from functools import partial
from multiprocessing import Pool
from pybedtools import BedTool


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enrichment between bed files.")

arg_parser.add_argument("region_file_1", help='bed file 1 (shuffled)')

arg_parser.add_argument("region_file_2", help='bed file 2 (not shuffled)')

arg_parser.add_argument("-i", "--iters", type=int, default=100,
                        help='number of simulation iterations; default=100')

arg_parser.add_argument("-s", "--species", type=str, default='hg19', choices=['hg19', 'hg38', 'mm10'],
                        help='species and assembly; default=hg19')

arg_parser.add_argument("-n", "--num_threads", type=int,
                        help='number of threads; default=SLURM_CPUS_PER_TASK or 1')

arg_parser.add_argument("--print_counts_to", type=str, default=None,
                        help="print expected counts to file")

arg_parser.add_argument("--elem_wise", action='store_true', default=False,
                        help='perform element-wise overlaps; default=False')

args = arg_parser.parse_args()

# save parameters
ANNOTATION_FILENAME = args.region_file_1
TEST_FILENAME = args.region_file_2
COUNT_FILENAME = args.print_counts_to
ITERATIONS = args.iters
SPECIES = args.species
ELEMENT = args.elem_wise

# calculate the number of threads
if args.num_threads:
    num_threads = args.num_threads
else:
    num_threads = int(os.getenv('SLURM_CPUS_PER_TASK', 1))


###
#   functions
###
def loadConstants(species):
    return {'hg19': ("/dors/capra_lab/users/bentonml/data/hg19_blacklist_gap.bed", "/dors/capra_lab/data_clean/dna/human/hg19/hg19_trim.chrom.sizes"),
            'hg38': ("/dors/capra_lab/users/bentonml/data/hg38_blacklist_gap.bed", "/dors/capra_lab/data_clean/dna/human/hg38/hg38_trim.chrom.sizes"),
            'mm10': ("/dors/capra_lab/users/bentonml/data/mm10_blacklist_gap.bed", "/dors/capra_lab/data_clean/dna/mouse/mm10/mm10_trim.chrom.sizes")
            }[species]


def calculateObserved(annotation, test, elementwise):
    obs_sum = 0

    if elementwise:
        obs_sum = annotation.intersect(test, u=True).count()
    else:
        obs_intersect = annotation.intersect(test, wo=True)
        for line in obs_intersect:
            obs_sum += int(line[-1])

    return obs_sum


def calculateExpected(annotation, test, elementwise, species, iters):
    BLACKLIST, CHROM_SZ = loadConstants(species)
    exp_sum = 0

    rand_file = annotation.shuffle(genome='hg19', excl=BLACKLIST, chrom=True, noOverlapping=True)

    if elementwise:
        exp_sum = rand_file.intersect(test, u=True).count()
    else:
        exp_intersect = rand_file.intersect(test, wo=True)
        for line in exp_intersect:
            exp_sum += int(line[-1])

    return exp_sum


def calculateEmpiricalP(obs, exp_sum_list):
    mu = np.mean(exp_sum_list)
    sigma = np.std(exp_sum_list)
    dist_from_mu = [exp - mu for exp in exp_sum_list]
    p_sum = sum(1 for exp_dist in dist_from_mu if abs(exp_dist) >= abs(obs - mu))

    fold_change = (obs + 1.0) / (mu + 1.0)
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
    obs_sum = calculateObserved(ANNOTATION_FILENAME, TEST_FILENAME, ELEMENT)

    # create pool and run simulations in parallel
    pool = Pool(num_threads)
    partial_calcExp = partial(calculateExpected, ANNOTATION_FILENAME, TEST_FILENAME, ELEMENT, SPECIES)
    exp_sum_list = pool.map(partial_calcExp, [i for i in range(ITERATIONS)])

    # wait for results to finish before calculating p-value
    pool.close()
    pool.join()

    # calculate empirical p value
    print(calculateEmpiricalP(obs_sum, exp_sum_list))

    if COUNT_FILENAME is not None:
        with open(COUNT_FILENAME, "w") as count_file:
            count_file.write('{}\n{}\n'.format(obs_sum, '\t'.join(map(str, exp_sum_list))))


if __name__ == "__main__":
    main(sys.argv[1:])

