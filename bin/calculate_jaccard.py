#!/bin/python
#
# marylaurenbenton | 2018
#
# this script will calculate the jaccard and relative jaccard between SORTED bed files
#
# relative_jaccard = jaccard / max_jaccard
# where: jaccard = |intersect| / ( |a| + |b| - |intersect| )
#        max_jaccard = min( |a|, |b| ) / max( |a|, |b| )
#
# pass in the names of 2 bed files & optional argument to specify number of significant digits
#

import sys
import argparse
import pybedtools


###
#  arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate Jaccard and relative Jaccard similarity between bed files.")

arg_parser.add_argument("bed_file_1", help='first BED file; should be sorted')
arg_parser.add_argument("bed_file_2", help='second BED file; should be sorted')
arg_parser.add_argument('-d', '--decimal', type=int, default=3, help='number of decimal places | default = 3')

args = arg_parser.parse_args()

# save parameters
A = args.bed_file_1
B = args.bed_file_2
DECIMAL = args.decimal


###
#  main
###
def main(argv):
    # create bedtools objects
    a = pybedtools.BedTool(A)
    b = pybedtools.BedTool(B)

    # calculate jaccard similarity
    result = a.jaccard(b)

    # calculate len of A and B
    a_self = a.jaccard(a)['intersection']
    b_self = b.jaccard(b)['intersection']

    # calculate the relative jaccard
    relative = result['jaccard'] / (min(a_self, b_self) / max(a_self, b_self))

    print('Jaccard: {}'.format(round(result['jaccard'], DECIMAL)))
    print('Relative Jaccard: {}'.format(round(relative, DECIMAL)))


if __name__ == "__main__":
    main(sys.argv[1:])

