#!/bin/python
#
# Mary Lauren Benton, 2017
#
# Calculate empirical p-value to compare enrichment between simulations
# Assumes comparison of region from the shuffled files; performs t-test
#

import sys
import argparse
import datetime
import numpy as np
from scipy import stats

###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Compare enrichment between count files.")

arg_parser.add_argument("count_file_1", help="count file for region_file_1")
arg_parser.add_argument("count_file_2", help="count file for region_file_2")

args = arg_parser.parse_args()

# save parameters
COUNTS_A_FILENAME = args.count_file_1
COUNTS_B_FILENAME = args.count_file_2


###
#   functions 
###
def read_counts(filename):
    with open(filename, 'r') as infile:
        obs = int(infile.readline().strip())

        # expects data in counts format: [ obs ]
        #                                [ exp ]
        #                                [ ... ]
        exp_list = []
        for exp in infile:
            exp_list.append(exp.strip())

    return obs, np.array(exp_list, dtype=int)


def calculate_fc_dist(obs, exp_list):
    return np.array([(obs + 1.0) / (exp + 1.0) for exp in exp_list])


###
#   main 
###
def main(argv):
    # print header
    print('{:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())[:20]))
    print('TestStat\tp-value')

    obs_a, exp_list_a = read_counts(COUNTS_A_FILENAME)
    obs_b, exp_list_b = read_counts(COUNTS_B_FILENAME)

    fc_list_a = calculate_fc_dist(obs_a, exp_list_a)
    fc_list_b = calculate_fc_dist(obs_b, exp_list_b)

    result = stats.ttest_ind(fc_list_a, fc_list_b, equal_var=False)
    print('{:.3f}\t{:.3f}'.format(result.statistic, result.pvalue))


if __name__ == "__main__":
    main(sys.argv[1:])

