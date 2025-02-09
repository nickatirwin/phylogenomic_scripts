# protein_prediction.py 
# 
# predict proteins from a genome using miniprot (https://github.com/lh3/miniprot?tab=readme-ov-file)
#
# usage: python protein_prediction.py [proteins] [genome] [genetic code] [threads]
#
# note: genome file name should start with taxonomy id separated by period - eg. 9606.genome.fasta
#
# requirements:
# mamba install bioconda::miniprot

import sys
import subprocess

proteins = sys.argv[1]
genome = sys.argv[2]
gc = sys.argv[3]
threads = sys.argv[4]

# run miniprot
subprocess.call('miniprot -Iut16 -t ' + threads + ' --trans ' + genome + ' ' + proteins + ' -T ' + gc + ' > '+genome.split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot',shell=True)

# parse output
mini = open(genome.split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot','r').readlines()
out = open(genome.split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot.fasta','w')
n = 0
p = 1
while n < len(mini):
    out.write('>'+genome.split('.')[0]+'.'+mini[n].split('\t')[0]+'_'+str(p)+'\n'+mini[n+1].strip().split('\t')[1]+'\n')
    n += 2
    p += 1
out.close()
