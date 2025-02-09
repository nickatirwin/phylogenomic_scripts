# FILTER ALIGNMENT BY PERCENT POSITIONS
#
# Filter an alignment for a given amount of missing data.
# Returns a filtered alignment file in fasta format. The proportion is the percent of data required to not be filtered out.
# Output is *.filt_[PROPORTION].fasta
#
# Usage: python filter_alignment_by_percent_positions.py [alignment.fasta] [proportion]
# Example: python filter_alignment_by_percent_positions.py myalignment.fasta 0.2 
#

import sys
from glob import glob

#calculate lengths and write output

prop = float(sys.argv[2])


for i in glob(sys.argv[1]):
    n = 0
    infile = open(i, 'r').read()
    seqs = infile.split('>')[1:]
    out = open(i + '.cov_' + str(prop), 'w')
    while n < len(seqs):
        seq = seqs[n].split('\n')
        total_length = 0
        present_data = 0
        for line in seq[1:]:
            total_length = total_length + len(line)
            present_data = present_data + len(line) - line.count('-') - line.count('X')
        if float(present_data)/total_length >= prop:
            out.write('>'+seq[0]+'\n')
            for parts in seq[1:-1]:
                out.write(parts + '\n')
        else:
            print(seq[0] + ' deleted')
        n = n + 1
    out.close()
