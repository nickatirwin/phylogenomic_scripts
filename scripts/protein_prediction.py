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
subprocess.call('miniprot -Iut16 -t ' + threads + ' --trans ' + genome + ' ' + proteins + ' -T ' + gc + ' > '+genome.split('/')[-1].split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot',shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

# parse output
mini = open(genome.split('/')[-1].split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot','r').readlines()

if len(mini) > 0:
    out = open(genome.split('/')[-1].split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot.fasta','w')
    n = 0
    p = 1
    while n < len(mini):
        if '##STA' in mini[n+1]:
            out.write('>'+genome.split('/')[-1].split('.')[0]+'.seq_'+str(p)+'\n'+mini[n+1].strip().split('\t')[1]+'\n')
            p += 1
        n += 2
    out.close()

    # cluster the resulting sequences at 100% identity
    subprocess.call('cd-hit -i '+genome.split('/')[-1].split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot.fasta'+' -o ' + genome.split('/')[-1].split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot.fasta.100c' + ' -c 1.0 -T '+threads,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    # remove temp files
    subprocess.call('mv '+genome.split('/')[-1].split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot.fasta.100c '+genome.split('/')[-1].split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot.fasta',shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    subprocess.call('rm '+genome.split('/')[-1].split('.')[0]+'.'+proteins.split('.')[0]+'.miniprot.fasta.100c.clstr',shell=True)


