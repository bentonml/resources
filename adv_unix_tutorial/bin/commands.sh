#!/bin/bash
# mary lauren benton | 2018
# advanced unix / bash tutorial

file='../data/E073_18_core_K27ac_mnemonics.bed'
gwas='../data/gwasCatalog_2018-04-01_hg38_n-200.bed'

###
#  the classics
###
cat  $file        # print the entire file

head $file        # print first 10 lines of file
tail $file        # print last 10 lines of file

less $file        # view contents of file (up/down)
less -U $file     # shows tabs

wc -l $file       # count the number of lines in file
cut -f 2-3 $file  # cut columns 2-3 from file

sort $file        # sort the file
uniq $file        # print uniq lines in file


###
#  combine commands | with pipes!
###
cut -f4 $file | sort -un  # cut the 4th column of the file, sort unique lines numerically

sort -k1,1 -k2,2n $file | cut -f1-3 | uniq  # sort bed file by chromosome, then by numeric
                                            # start and find uniq coordinates

# detour!

# we could also do that by:
cut -f1-3 $file | sort -u -k1,1 -k2,2n

# and test our theory by:
sort -k1,1 -k2,2n $file | cut -f1-3 | uniq > tmp.cmd1
cut -f1-3 $file | sort -u -k1,1 -k2,2n     > tmp.cmd2
diff -q tmp.cmd1 tmp.cmd2  # -q supresses output besides 'differ' or 'don't differ'

# end detour.


cut -f1 $file | sort -u | wc -l  # count the number of unique entries in a column


###
#  searching for strings with grep
###
grep rs139770682 $gwas     # select lines with a match to this rsID
grep -v rs139770682 $gwas  # select lines that don't match
grep -i RS139770682 $gwas  # ignore case

grep -c European $gwas       # count matching lines
grep European $gwas | wc -l  # see? same thing!

echo -e 'Breast cancer\ngut microbiota\nheight' > test_patterns
grep -f test_patterns $gwas    # take patterns/regex from a file
grep -cf test_patterns $gwas   # you can combine options
grep -cif test_patterns $gwas  # you can combine options

grep -E '[^NR,]\sEuropean(;European)*?\s' $gwas  # interpret pattern as regex


###
#  regex: let's make it more advanced
###
grep -E '[^NR,]\sEuropean(;European)*?\s' $gwas  # created example with regex101.com

# . stands for any character
# * is 0 or more; + is one or more
# ? makes 'lazy' quantifiers
# () denote capture groups; [] denote character classes
# \s for whitespace characters; \d for digits


###
#  sed & awk — languages of their own...
###
sed '1p' $gwas     # print first line (twice b/c default printing)
sed -n '1p' $gwas  # supress default printing, print 1st line

sed '1d' $gwas     # print file without first line

sed 's/ \+/\t/' $gwas         # convert spaces to tabs
sed -i.bak 's/ \+/\t/' $gwas  # convert spaces to tabs, in place, create backup

sed 's/;/,/g' $gwas | head    # substitute ; for , - everywhere (all matches)
sed 's/European/AnotherEuropean/2p' $gwas | head -20  # make substitution on 2nd match only


awk '{print}' $gwas  # this is a basic awk script!

awk -F'\t' '{print $7}' $gwas  # use only tabs, not tabs AND spaces

awk '{$1=="chr13"; print}' $gwas  # print lines with $1 == chr13 only
awk '$1 ~ /chr13/' $gwas          # an alternate approach

awk '$1 ~ /chr13/ {print $1,"with start (", $2, ") and end (", $3, ")"}' $gwas  # adding text to printed output

# internal variables and extra blocks
awk -v OFS='\t' '$1 ~ /chr13/ {print $1, $2, $3, $3-$2}' $gwas      # prints region and length, tab-separated
awk 'BEGIN{OFS="\t"} $1 ~ /chr13/ {print $1, $2, $3, $3-$2}' $gwas  # alternatively
awk 'BEGIN{print "number of bp in this BED file:"} {sum += $3 - $2} END{print sum}' $gwas  # sample use of BEGIN and END blocks

# what if we want to use a variable?
match_this="chr13"
awk -v var="$match_this" '$1 ~ var {print}' $gwas  # note: we removed the '//' around the pattern

# boolean combinations are cool
awk '$1 ~ /chr13/ || $4 ~ /4E-9/ {print $1, $2, $3, $4}' $gwas  # on chromosome 13 OR p = 4E-9
awk -v FS='\t' '$1 ~ /chr13/ && $6 > 5000' $gwas  # on chromosome 13 AND with sample size > 5000


###
#  scripting in bash
###
# for details & examples, check out the standalone scripts in bin/


###
#  bonus! the bash_profile
###
ls -a  # to see hidden files like .bashrc and .bash_profile

# note: the bash_profile is for login shells (what you ssh into on ACCRE)
# the bashrc is for interactive non-login shells
# fun fact: OS X terminals always run login shells by default

