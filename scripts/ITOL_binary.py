# ITOL presence-absence plot
#
# usage: python ITOL_binary.py '*fasta'
#
# where fasta = geneName.fasta and headers = >taxid.proteinID

import sys
import subprocess
import random
from glob import glob

sep = ','

# 1. download the annotation file
subprocess.call('wget https://itol.embl.de/help/dataset_binary_template.txt',shell = True)
subprocess.call('mv dataset_binary_template.txt ITOL_binary_annotation.txt', shell = True)

# 2. change the optional values to fit the dataset
genes = [fname.split('.')[0] for fname in glob(sys.argv[1])]

# random colour selection
colours = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(0,len(genes))]

# change the title
subprocess.call("sed -i 's/DATASET_LABEL,label1/DATASET_LABEL,protein_prescence/g' ITOL_binary_annotation.txt",shell = True)

# redefine the legend
subprocess.call("sed -i 's/#LEGEND_TITLE,Dataset legend/LEGEND_TITLE,Dataset legend/g' ITOL_binary_annotation.txt", shell = True)
subprocess.call("sed -i 's/#LEGEND_SHAPES,1,2,3/LEGEND_SHAPES,"+sep.join(['1' for n in range(0,len(genes))])+"/g' ITOL_binary_annotation.txt", shell = True)
subprocess.call("sed -i 's/#LEGEND_LABELS,value1,value2,value3/LEGEND_LABELS,"+sep.join([g for g in genes])+"/g' ITOL_binary_annotation.txt", shell = True)
subprocess.call("sed -i 's/#LEGEND_COLORS,#ff0000,#00ff00,#0000ff/LEGEND_COLORS,"+sep.join([c for c in colours])+"/g' ITOL_binary_annotation.txt", shell = True)

# redefine the fields
subprocess.call("sed -i 's/FIELD_SHAPES,1/FIELD_SHAPES,"+sep.join(['1' for n in range(0,len(genes))])+"/g' ITOL_binary_annotation.txt", shell = True)
subprocess.call("sed -i 's/FIELD_LABELS,f1/FIELD_LABELS,"+sep.join([g for g in genes])+"/g' ITOL_binary_annotation.txt", shell = True)
subprocess.call("sed -i 's/#FIELD_COLORS,#ff0000/FIELD_COLORS,"+sep.join([c for c in colours])+"/g' ITOL_binary_annotation.txt", shell = True)

# collect data for each gene
# record all taxa
t = []
for fname in glob(sys.argv[1]):
    f = open(fname,'r').readlines()
    for line in f:
        if line.startswith('>'):
            t.append(line.split('.')[0].strip('>'))
t = list(set(t))

taxa_d = {}
for taxa in t:
    taxa_d[taxa] = ['-1' for n in range(0,len(genes))]

# record data
n = 0
for fname in glob(sys.argv[1]):
    f = open(fname,'r').readlines()
    for line in f:
        if line.startswith('>'):
            taxa_d[line.split('.')[0].strip('>')][n] = '1'
    n += 1

# output the results
out = open('ITOL_binary_annotation.txt','a')
for t in taxa_d:
    out.write(t+','+sep.join(taxa_d[t])+'\n')
out.close()


