# FILTER SEQUENCE LENGTHS
#
# Filter a fasta file for a given sequence length
# Returns a filtered fasta file.
# Output is *.filt_[LENGTH].fasta
#
# Usage: python filter_sequence_lengths.py [fasta] [min/max] [min_seq_length]
# Example: python filter_sequence_lengths.py '*.trimal' max 100

import sys
from glob import glob

#calculate lengths and write output

fasta = sys.argv[1]
m = str(sys.argv[2]).lower()
mlen = int(sys.argv[3])

seq_d = {}
f = open(fasta,'r').read().split('>')[1:]
for l in f:
    seq_d[l.split('\n')[0]] = [l.split('\n',1)[1].replace('\n',''),len(l.split('\n',1)[1].replace('\n',''))]
out = open(fasta.rsplit('.',1)[0]+'.filt_'+str(mlen)+'.fasta','w')
if 'min' in m:
    for s in seq_d:
        if seq_d[s][1] >= mlen:
            out.write('>'+s+'\n'+seq_d[s][0]+'\n')
elif 'max' in m:
    for s in seq_d:
        if seq_d[s][1] <= mlen:
            out.write('>'+s+'\n'+seq_d[s][0]+'\n')
out.close()


