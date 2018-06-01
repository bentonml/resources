#!/bin/python
#
# Mary Lauren Benton, 2018
#

import re
import datetime

###
#  FUNCTIONS
###

# returns list with text from all capture groups
def match_to_string(m):
    return ' '.join([g.lower() for g in m.groups() if g is not None])

def print_to_file(codes):
    with open('./file.out', 'w') as outfile:
        outfile.write('###\n'
                      '#  author: [ name ] \n'
                      '#  date:   {}\n'
                      '#  output: file.out\n'
                      '###\n'
                      .format(str(datetime.datetime.now())[:10]))
        for c in sorted(codes): outfile.write(c + '\n')

###
#  MAIN
###
pattern = 'regex here'
regex   = re.compile(pattern, flags=re.IGNORECASE) # if case sensitive remove flag

results = []

with open('file.dat', 'r') as data_file:
    for line in data_file:
        m = regex.search(line)  # finds all occurences
        if m is not None:
            # can take any action with line/match here
            results.add(match_to_string(m))

print_to_file(results)

