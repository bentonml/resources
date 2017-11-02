#!/bin/python
#
# Mary Lauren Benton, 2017
#
# Depends on:
#               BEDtools v2.23.0-20
#               /dors/capra_lab/data/dna/[species]_trim.chrom.sizes
#               /dors/capra_lab/bentonml/data/blacklist_gap_[species].bed
#

import sys
import os
import subprocess
import argparse
import datetime
import numpy as np
from functools import partial
from multiprocessing import Pool 


def loadConstants():
    return "/dors/capra_lab/bentonml/pavlicev_atacseq_enrichment/data/blacklist_gap_mm10.bed", "/dors/capra_lab/data/dna/mm10_trim.chrom.sizes"


def calculateObserved(annotation, test, elementwise):
    obs_sum = 0

    if elementwise:
        obs_intersect = subprocess.check_output(['intersectBed', '-u', '-a', annotation, '-b', test])
        obs_sum = len(obs_intersect.splitlines())
    else:
        obs_intersect = subprocess.check_output(['intersectBed', '-wo', '-a', annotation, '-b', test])

        for line in obs_intersect.splitlines():
            obs_sum += int(line.split('\t')[-1])

    return obs_sum


def calculateExpected(annotation, test, elementwise, iters):
    BLACKLIST, CHROM_SZ = loadConstants()
    exp_sum = 0

    rand_file = subprocess.Popen(['shuffleBed', '-excl', BLACKLIST, '-i', annotation, '-g', CHROM_SZ, '-chrom', '-noOverlapping'], stdout=subprocess.PIPE)

    if elementwise:
        exp_intersect = subprocess.check_output(['intersectBed', '-u', '-a', 'stdin', '-b', test], stdin=rand_file.stdout)
        exp_sum = len(exp_intersect.splitlines())
    else:
        exp_intersect = subprocess.check_output(['intersectBed', '-wo', '-a', 'stdin', '-b', test], stdin=rand_file.stdout)

        for line in exp_intersect.splitlines():
            exp_sum += int(line.split('\t')[-1])

    return exp_sum


def calculateEmpiricalP(obs, exp_sum_list):
    mu = np.mean(exp_sum_list)
    dist_from_mu = [exp - mu for exp in exp_sum_list]
    p_sum = sum(1 for exp_dist in dist_from_mu if abs(exp_dist) >= abs(obs - mu))

    fold_change = (obs + 1.0) / (mu + 1.0)
    p_val = (p_sum + 1.0) / (len(exp_sum_list) + 1.0) 

    return "%d\t%.3f\t%.3f\t%.3f" % (obs, exp, fold_change, p_val)


def usage():
    print "usage: python calculate_enrichment.py [annotation file] [test file] [iterations] [element-wise T/F]"


def main(argv):
    # check only for correct number of parameters
    if len(argv) < 3 or len(argv) > 5:
        usage()
        sys.exit(1)
    
    # save parameters
    ANNOTATION_FILENAME= sys.argv[1]
    TEST_FILENAME = sys.argv[2]
    ITERATIONS = int(sys.argv[3])

    if sys.argv[4] == "T":
        ELEMENT = True
    else:
        ELEMENT = False

    # calculate the number of threads
    if len(argv) == 5:
        num_threads = int(sys.argv[5])
    else:
        num_threads = int(os.environ['SLURM_CPUS_PER_TASK'])

    # print header
    print('{:s} {}'.format(' '.join(sys.argv), str(datetime.datetime.now())[:20]))
    print('Observed\tExpected\tFoldChange\tp-value')
    
    # run initial intersection and save
    obs_sum = calculateObserved(ANNOTATION_FILENAME, TEST_FILENAME, ELEMENT)

    # create pool and run simulations in parallel
    pool = Pool(num_threads)
    partial_calcExp = partial(calculateExpected, ANNOTATION_FILENAME, TEST_FILENAME, ELEMENT)
    exp_sum_list = pool.map(partial_calcExp, (i for i in range(ITERATIONS)))
    
    # wait for results to finish before calculating p-value
    pool.close()
    pool.join()

    # calculate empirical p value 
    print calculateEmpiricalP(obs_sum, exp_sum_list)


if __name__ == "__main__":
    main(sys.argv[1:])

