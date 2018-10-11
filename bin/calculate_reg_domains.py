###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2018.06.25
#   updated   | 2018.06.26
#
#   this script will calculate a 'basal plus extension' regulatory domain using
#   the definitions from GREAT and save the resulting bed file
###


import sys
import pybedtools
import argparse
import pandas as pd
import numpy  as np


###
#  variables
###
arg_parser = argparse.ArgumentParser(description="Create a BED file with the regulatory domains for a set of genes.")

arg_parser.add_argument('gene_file',   help='BED file of genes; should be sorted')
arg_parser.add_argument('-a', '--algorithm',  type=str, default='great', choices=['great', 'window'], help='regulatory domain definition; default=great')
arg_parser.add_argument('-s', '--species',    type=str, default='hg19', choices=['hg19', 'hg38'], help='species and assembly; default=hg19')
arg_parser.add_argument('-u', '--upstream',   type=int, default=5, help='basal upstream extension in kb; default=5')
arg_parser.add_argument('-d', '--downstream', type=int, default=1, help='basal downstream extension in kb; default=1')
arg_parser.add_argument('-e', '--extension',  type=int, default=1000, help='maximum extension in kb; default=1000')
arg_parser.add_argument('-o', '--outfile',    type=str, default='result.bed', help='output file name; default=result.bed')


args = arg_parser.parse_args()

# save parameters
UP_EXTENSION = args.upstream   * 1000
DN_EXTENSION = args.downstream * 1000
MAX_EXTENSION = args.extension * 1000
SPECIES = args.species
GENE_FILE = args.gene_file
ALGORITHM = args.algorithm
OUTFILE = args.outfile


###
#  functions
###
# modified from https://stackoverflow.com/questions/323750/how-to-access-previous-next-element-while-for-looping
def neighbors(iterable):
    itr = iter(iterable)
    prev_item = None
    curr_item = next(itr)
    for next_item in itr:
        yield (prev_item, curr_item, next_item)
        prev_item, curr_item = curr_item, next_item
    yield (prev_item, curr_item, None)


def basal_plus_extension(genes, species, up_extension, dn_extension, max_extension):
    chrom_sizes = pybedtools.helpers.chromsizes(species)
    expanded_genes = genes.slop(genome=species, l=(up_extension), r=(dn_extension), s=True)
    result = ''

    for p, c, n in neighbors(expanded_genes):
        chr_size = chrom_sizes[c.chrom][1]
        curr_tss = c.start + up_extension if c.strand == '+' else c.start + dn_extension

        tmp_start = max(0, curr_tss - max_extension)
        basal_start = c.start
        tmp_start = min(basal_start, tmp_start)

        if (p is not None) and (p.chrom == c.chrom):
            prev_tss = p.start + up_extension if p.strand == '+' else p.start + dn_extension
            prev_end = prev_tss + dn_extension if p.strand == '+' else prev_tss + up_extension
            tmp_start = min(basal_start, max(prev_end, tmp_start))

        tmp_end = min(chr_size, curr_tss + max_extension);
        basal_end = c.end;
        tmp_end = max(basal_end, tmp_end)

        if (n is not None) and (c.chrom == n.chrom):
            next_tss = n.start + up_extension if n.strand == '+' else n.start + dn_extension
            next_start = next_tss - up_extension if n.strand == '+' else next_tss - dn_extension
            tmp_end = max(basal_end, min(next_start, tmp_end))

        result += '\n' + c.chrom + ' ' + str(tmp_start) + ' ' + str(tmp_end) + ' ' + c.name + ' ' + c.strand

    return result


def gene_window(genes, species, extension):
    return genes.flank(genome=species, b=extension)


###
#  main
###
def main(argv):
    # save variables
    genes = pybedtools.BedTool(GENE_FILE)

    if ALGORITHM == 'great':
        pybedtools.BedTool(basal_plus_extension(genes, SPECIES, UP_EXTENSION, DN_EXTENSION, MAX_EXTENSION),
                           from_string=True).saveas(OUTFILE)
    elif ALGORITHM == 'window':
        gene_window(genes, SPECIES, MAX_EXTENSION).saveas(OUTFILE)
    else:
        print('Invalid algorithm option.')


if __name__ == "__main__":
    main(sys.argv[1:])

