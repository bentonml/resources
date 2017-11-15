#!/bin/python
#
# Mary Lauren Benton, 2017
#
# Calculate empirical p-value to compare enrichment between simulations
# Assumes comparison of region from the shuffled files
#

import sys
import argparse
import datetime
import numpy as np

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

        # expects data in counts format: [ obs ] [ exp \t exp \t ... ]
        exp_list = infile.readline().strip().split('\t')

    return obs, np.array(exp_list, dtype=int)


def calculate_fc_distribution(obs, exp_list):
    return np.array([(obs + 1.0) / (exp + 1.0) for exp in exp_list])


def calculate_empirical_p(obs_fc_a, fc_list_a, obs_fc_b, fc_list_b):
    delta_final = obs_fc_a - obs_fc_b
    delta_dist  = [fc_list_a[i] - fc_list_b[i] for i in range(len(fc_list_a))]

    # sum number of differences >= to observed difference
    p_sum = sum(1 for delta in delta_dist if abs(delta) >= abs(delta_final))

    return (p_sum + 1.0) / (len(delta_dist) + 1.0)


###
#   main 
###
def main(argv):
    # print header
    print('{:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())[:20]))
    print('p-value')

    obs_a, exp_list_a = read_counts(COUNTS_A_FILENAME)
    obs_b, exp_list_b = read_counts(COUNTS_B_FILENAME)

    fc_list_a = calculate_fc_distribution(obs_a, exp_list_a)
    fc_list_b = calculate_fc_distribution(obs_b, exp_list_b)

    final_obs_a = (obs_a + 1.0) / (np.mean(exp_list_a) + 1.0)
    final_obs_b = (obs_b + 1.0) / (np.mean(exp_list_b) + 1.0)

    p_val = calculate_empirical_p(final_obs_a, fc_list_a, final_obs_b, fc_list_b)

    # print result
    print('{:.3f}'.format(p_val))


if __name__ == "__main__":
    main(sys.argv[1:])

