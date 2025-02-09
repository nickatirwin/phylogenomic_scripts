# ITOL heatmap plot
#
# usage: python ITOL_heatmap.py '*fasta'
#
# where fasta = geneName.fasta and headers = >taxid.proteinID

import sys
import subprocess
import random
from glob import glob
from collections import Counter

sep = ' '

# 1. download the annotation file
subprocess.call('wget https://itol.embl.de/help/dataset_heatmap_template.txt',shell = True)
subprocess.call('mv dataset_heatmap_template.txt ITOL_heatmap_annotation.txt', shell = True)

# 2. change the optional values to fit the dataset
genes = [fname.split('.')[0] for fname in glob(sys.argv[1])]

# colour gradient selection
min = '#deebf7'
max = '#3182bd'
colours = [min,max]
# change the title
subprocess.call("sed -i 's/DATASET_LABEL,label1/DATASET_LABEL,paralog_abundance/g' ITOL_heatmap_annotation.txt",shell = True)

# redefine the fields
subprocess.call("sed -i 's/#COLOR_MIN #ff0000/COLOR_MIN "+min+"/g' ITOL_heatmap_annotation.txt", shell = True)
subprocess.call("sed -i 's/#COLOR_MAX #0000ff/COLOR_MIN "+max+"/g' ITOL_heatmap_annotation.txt", shell = True)
subprocess.call("sed -i 's/FIELD_LABELS f1 f2 f3 f4 f5 f6/FIELD_LABELS " + sep.join([g for g in genes]) + "/g' ITOL_heatmap_annotation.txt", shell = True)

# collect data for each gene
# record all taxa
t = []
for fname in glob(sys.argv[1]):
    f = open(fname,'r').read().split('>')[1:]
    for line in f:
        t.append(line.split('.')[0].strip())
t = list(set(t))

taxa_d = {}
for taxa in t:
    taxa_d[taxa.strip()] = ['0' for n in range(0,len(genes))]

# record data
n = 0
for fname in glob(sys.argv[1]):
    f = open(fname,'r').read().split('>')[1:]
    counts = Counter([t.split('.')[0].strip() for t in f])
    for t in taxa_d:
        try:
            taxa_d[t][n] = str(counts[t])
        except:
            pass
    n += 1

# output the results
out = open('ITOL_heatmap_annotation.txt','a')
for t in taxa_d:
    out.write(t.strip()+' '+sep.join(taxa_d[t])+'\n')
out.close()
