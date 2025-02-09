# Remove illegal characters from protein fasta
#
# Usage: python remove_illegal.py [fasta]
# Eg. python remove_illegal.py proteinA.fasta

import sys

########################################
# define illegal amino acid identifiers
illegal_aa = ['U','O','J','Z','B']
# replace illegal amino acids with:
replace = 'X'
########################################

fasta = open(sys.argv[1],'r').read().split('>')[1:]

out = open(sys.argv[1]+'.clean','w')

for seq in fasta:
    seq = seq.split('\n')
    header = '>'+seq[0].strip()
    sequence = ''
    for part in seq[1:]:
        sequence = sequence + part.strip()
    for c in illegal_aa:
        if c in sequence:
            sequence = sequence.replace(c,replace)
    out.write(header+'\n'+sequence+'\n')
    
out.close()
