#!/bin/python
###
#   name    | mary lauren benton
#   created | 2017
#   updated | 2021
#
#   depends on:
#       BEDtools v2.23.0-20 via pybedtools
###

import os
import sys, traceback
import argparse
import numpy as np
from pybedtools import BedTool
from pybedtools.helpers import BEDToolsError, cleanup, get_tempdir, set_tempdir


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Generate length-matched random bed file.")

arg_parser.add_argument("input_bed", help='bed file to randomize')


arg_parser.add_argument("-s", "--species", type=str, default='hg19', choices=['hg19', 'hg38'],
                        help='species and assembly; default=hg19')

arg_parser.add_argument("-b", "--blacklist", type=str, default=None,
                        help='custom blacklist file; default=None')

arg_parser.add_argument("-o", "--outfile", type=str, default='random.bed',
                        help="name of randomized bed file; default=random.bed")


args = arg_parser.parse_args()

# save parameters
INPUT_FILENAME = args.input_bed
SPECIES = args.species
BLACKLIST = args.blacklist
OUTPUT_FILENAME = args.output


# if running on slurm, set tmp to runtime dir
set_tempdir(os.getenv('ACCRE_RUNTIME_DIR', get_tempdir()))


###
#   main
###
def main(argv):

    # generate random bed file based on input bed file
    try:
        rand_file = INPUT_FILENAME.shuffle(genome=SPECIES, excl=BLACKLIST, 
                                           chrom=True, noOverlapping=True)
    except BEDToolsError:
        print('ERROR: Shuffling produced BEDToolsError.')
        cleanup()
        exit(1)

    # write ouput file
    rand_file.moveto(OUTPUT_FILENAME)


    # clean up any pybedtools tmp files
    cleanup()


if __name__ == "__main__":
    main(sys.argv[1:])

