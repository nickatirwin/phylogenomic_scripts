# single_line_fasta.py
#
# Convert a multiline fasta into a single line fasta.
# Outputs *.fasta.sl
#
# Usage: python single_line_fasta.py [multiline.fasta]
# Example: python single_line_fasta.py multiline.fasta
#

import sys

n = 0

f = open(sys.argv[1], 'r').read().split('>')[1:]
output = open(sys.argv[1] + '.sl', 'w')

for seq in f:
    output.write('>'+seq.split('\n')[0]+'\n'+seq.split('\n',1)[1].replace('\n','')+'\n')
output.close()
