#!/usr/env/bin python
#
# Mary Lauren Benton, 2016
#
# This script will pull selected columns from an input file and generate a BED
# file for downstream analysis.

import sys
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='Parses input file into BED format', epilog="// mary lauren benton")
    reqgrp = parser.add_argument_group('required arguments')
    optgrp = parser.add_argument_group('other positional arguments')

    reqgrp.add_argument("-i", "--input", help="name of input file", nargs='?', dest="infile",  type=argparse.FileType('r'), required=True)
    parser.add_argument("-o", "--output", help="name of output file (default: stdout)", nargs='?', dest="outfile", type=argparse.FileType('w'), default=sys.stdout)

    reqgrp.add_argument("chrom", type=int, help="idx of column with chromosome name", default=-1)
    reqgrp.add_argument("start", type=int, help="idx of column with start coordinate", default=-1)
    reqgrp.add_argument("end",   type=int, help="idx of column with end coordinate", default=-1)

    optgrp.add_argument("name",   type=str, nargs='?', help="short name/label for file ('.' if N/A)", default='.')
    optgrp.add_argument("score",  type=int, nargs='?', help="idx of column with score (-1 if N/A)", default=-1)
    optgrp.add_argument("strand", type=int, nargs='?', help="idx of column with strand (+/-, or -1 if N/A)", default=-1)
    # better way to specify empty columns?

    return parser.parse_args()

def checkBounds(args, numCol):
    if not (0 < args.chrom <= numCol):
        print "Invalid index value (chrom): %d" % (args.chrom)
        sys.exit(1)
    elif not (0 < args.start <= numCol):
        print "Invalid index value (start): %d" % (args.start)
        sys.exit(1)
    elif not (0 < args.end <= numCol):
        print "Invalid index value (end): %d" % (args.end)
        sys.exit(1)
    elif args.score < -1 or args.score == 0 or args.score > numCol:
        print "Invalid index value (score): %d" % (args.score)
        sys.exit(1)
    elif args.strand < -1 or args.score == 0 or args.strand > numCol:
        print "Invalid index value (strand): %d" % (args.strand)
        sys.exit(1)

def idx(i): # assumes numbering from 1
    return i - 1

def generateBED(args, n):
    if args.score == -1 and args.strand == -1:
        args.outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(n[idx(args.chrom)], n[idx(args.start)], n[idx(args.end)], args.name, '.', '.'))
    elif args.score == -1:
        args.outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(n[idx(args.chrom)], n[idx(args.start)], n[idx(args.end)], args.name, '.', n[idx(args.strand)]))
    elif args.strand == -1:
        args.outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(n[idx(args.chrom)], n[idx(args.start)], n[idx(args.end)], args.name, n[idx(args.score)], '.'))
    else:
        args.outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(n[idx(args.chrom)], n[idx(args.start)], n[idx(args.end)], args.name, n[idx(args.score)], n[idx(args.strand)]))

def main():
    args = parseArguments()

    for line in args.infile:
        n = line.strip().split('\t')
        checkBounds(args, len(n)) # boundary checking for column values (assumes numbered from 1)
        generateBED(args, n)

if __name__ == "__main__":
    main()

