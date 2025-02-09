# concatenation.py
'''
Concatenate alignments together for phylogenomic analyses.
Sequence headers must start with the species identifier followed by
a period (e.g. >Homo_sapiens.PROTEINID).

Output is a concatenated alignment file (concatenation.DATE.fasta) and
an overview of gene/site presence per species (.species_stat.tab file)
and the percentage of species representation for each gene (.gene_stat.tab).

Usage: python concatenation.py '*fasta.aln'
'''

from glob import glob
import sys
from datetime import datetime
from collections import Counter

# record all species
species = []
for fname in glob(sys.argv[1]):
    fasta = open(fname,'r').read().split('>')[1:]
    species += [s.split('.')[0] for s in fasta]

species_d = {}
counts_d = {}
aln_records = {}

for sp in list(set(species)):
    species_d[sp] = ''
    counts_d[sp] = 0

# add sequence information
for aln in glob(sys.argv[1]):
    aln_d = {}
    aln_file = open(aln,'r').read().split('>')[1:]
    length = len(aln_file[0].split('\n',1)[1].replace('\n',''))
    aln_records[aln.split('/')[-1].split('.')[0]] = 0
    for seq in aln_file:
        aln_d[seq.split('.')[0]] = seq.split('\n',1)[1].replace('\n','')
        aln_records[aln.split('/')[-1].split('.')[0]] += 1
    for sp in species_d:
        try:
            species_d[sp] += aln_d[sp]
            counts_d[sp] += 1
        except:
            species_d[sp] += ('-'*length)

# output concatenated alignment
t = []
for sp in species_d:
    t.append(sp)
taxa_counts = dict(Counter(t))
values = set(list(taxa_counts.values()))
if len(values) > 1:
    print('\nError: there can only be one sequence per species in each alignment\n')
else:
    out = open('concatenation.'+str(datetime.now()).split(' ')[0].replace('-','_')+'.fasta','w')
    for sp in species_d:
        out.write('>'+sp+'\n'+species_d[sp]+'\n')
    out.close()

# output species presence/absence statistics
total_genes = len(glob(sys.argv[1]))
total_length = len(species_d[sp])

out = open('concatenation.'+str(datetime.now()).split(' ')[0].replace('-','_')+'.species_stats.tab','w')
out.write('species\t%genes\t%sites\n')
for sp in species_d:
    out.write(sp+'\t'+str(round(100*(float(counts_d[sp])/total_genes),2))+'\t'+str(round(100*((total_length - float(species_d[sp].count('-')))/total_length),2))+'\n')
out.close()

# output the gene representation data
total_sps = len(list(set(species)))
out = open('concatenation.'+str(datetime.now()).split(' ')[0].replace('-','_')+'.gene_stats.tab','w')
out.write('gene\t%species\n')
for gene in list(reversed(sorted(aln_records, key=aln_records.get))):
    out.write(gene+'\t'+str(round(100*(aln_records[gene]/float(total_sps)),2))+'%\n')
out.close()
    
