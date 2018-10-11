#!/bin/bash
# mary lauren benton | 2018
#
# proposed solutions for problems 1-5

###
#  1
###

# see spaces instead of tabs
less -U ../data/gwasCatalog_2018-04-01_hg38_n-200.bed

# convert 1 or more spaces to a tab (1st match per line)
sed -i.bak's/ \+/\t/' ../data/gwasCatalog_2018-04-01_hg38_n-200.bed


###
#  2
###

# search for rs139770682
grep rs139770682 ../data/gwasCatalog_2018-04-01_hg38_n-200.bed


###
#  3
###

# get european or african descent
egrep '(European|Afr)' ../data/gwasCatalog_2018-04-01_hg38_n-200.bed


###
#  4
###

# without header line
awk -v OFS='\t' '$4 ~ "Enh(A.|Wk)" {print $1,$2,$3,$4,".","."}' ../data/E073_18_core_K27ac_mnemonics.bed

# with header
awk 'BEGIN{OFS="\t"; print "#chr","start","end","score","strand"} $4 ~ "Enh(A.|Wk)" {print $1,$2,$3,$4,".","."}' ../data/E073_18_core_K27ac_mnemonics.bed


###
#  5
###

# use regex to pull out relevant lines, awk to rearrange, sed to remove extra characters
zcat ../data/Homo_sapiens.GRCh37.75.gtf.gz | awk 'BEGIN{OFS="\t"} $3 ~ /gene/ && $1 ~ /^[0-9XY]/ {print "chr"$1,$4-1,$5,$10,$16,$7}' | sed 's/[";]//g'

# alternately, tr will delete the extra " and ; characters
zcat ../data/Homo_sapiens.GRCh37.75.gtf.gz | awk 'BEGIN{OFS="\t"} $3 ~ /gene/ && $1 ~ /^[0-9XY]/ {print "chr"$1,$4-1,$5,$10,$16,$7}' | tr -d '";'

