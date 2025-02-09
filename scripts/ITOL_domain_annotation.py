###  MAKE ITOL DOMAIN ANNOTATION FROM HMMSCAN
###
###  15-07-2020
#
# Overview: Annotate a tree with domain/protein structures. Make an ITOL domain 
# annotation file from the output of HMMscan. Place your tree in ITOL and drag 
# and drop the resulting annotation file onto it. 
#
# Usage: python ITOL_domain_annotation.py [fasta file] [parsed hmmscan file]
# Example: python ITOL_domain_annotation.py proteinA.fasta ProteinA.fasta.hmmscan.parsed
#
# NOTE: If you use this with an IQtree, you may have to remove quotes from your iqtree if you had special
#       characters in your headersi (eg. @). Do this by: sed -i "s/'//g" proteinA.treefile 
#
# The parsed HMM scan file can be obtained by running the following commands:
#
# hmmscan -E 1e-5 --incE 1e-5 --domE 1e-5 --cpu 5 -o ProteinA.fasta.hmmscan Pfam-A.hmm ProteinA.fasta
# python /Data/nick/scripts/parse_HMMscan_alignments.py ProteinA.all.fasta.hmmscan

"""
EXAMPLE PIPELINE:

# download the Pfam database
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

# unzip the database
gunzip Pfam-A.hmm.gz

# compress the database for use with the HMMscan function
hmmpress Pfam-A.hmm

# annotate your proteins
hmmscan -E 1e-5 --incE 1e-5 --domE 1e-5 --cpu 5 -o ProteinA.all.fasta.hmmscan Pfam-A.hmm ProteinA.all.fasta

# parse the output file
python parse_HMMscan_alignments.py ProteinA.all.fasta.hmmscan

# make your annotation file (specify iqtree if you made your tree in iqtree)
python ITOL_domain_annotation.py ProteinA.all.fasta ProteinA.all.fasta.hmmscan.parsed iqtree

# remove apostrophes (only necessary if you have special characters in the headers and when using iqtree because it replaces special characters)
sed -i "s/'//g" ProteinA.all.fasta.treefile

# Upload your tree to ITOL and annotate with your ITOL annotation file.

"""

import sys
import subprocess
import random

# 1. Load in fasta file

fasta = open(sys.argv[1],'r').read().split('>')[1:]

# record proteins and their lengths

proteins = {}
for line in fasta:
    seq = line.split('\n')[0].strip('>').strip()
    length = str(len(line.split('\n',1)[1].strip().replace('\n','')))
    proteins[seq] = [length]

# 2. Load in HMMScan file and get domain information

all_domains = []
HMMscan = open(sys.argv[2],'r').readlines()[1:]
for line in HMMscan:
    seq = line.split('\t')[2].rsplit('_',1)[0]
    try:
        proteins[seq].append(line.split('\t')[11]) #start
        proteins[seq].append(line.split('\t')[12]) #stop
        proteins[seq].append(line.split('\t')[1]) #domain
        all_domains.append(line.split('\t')[1])
    except: # the sequence in the HMMScan file was not in the fasta file
        pass

# 2. Import ITOL domain annotation template file

subprocess.call('wget https://itol.embl.de/help/dataset_protein_domains_template.txt',shell = True)
subprocess.call('mv dataset_protein_domains_template.txt ' + sys.argv[1] + '.ITOL.txt',shell = True)

ITOL = open(sys.argv[1]+'.ITOL.txt','a')

# 4. Add annotations
# define colours randomly

all_domains = list(set(all_domains))

colour = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(len(all_domains))]

colours = {}
n = 0
for d in all_domains:
    colours[d] = colour[n]
    n += 1

# add annotations to ITOL annotation file
annotation = ''
for seq in list(proteins.keys()):
    if len(proteins[seq]) > 1:
        annotation = annotation + seq + ',' + proteins[seq][0] + ','
    number_of_domains = (len(proteins[seq])-1)/3
    n = 0
    i = 0
    while i < number_of_domains:
        start = proteins[seq][n+1]
        stop = proteins[seq][n+2]
        domain = proteins[seq][n+3]
        colour = colours[domain]
        domain_annot = 'RE|'+start+'|'+stop+'|'+colour+'|'+domain+','
        annotation = annotation + domain_annot
        n += 3
        i += 1
    annotation = annotation.strip(',')+'\n'

while '\n\n' in annotation:
    annotation = annotation.replace('\n\n','\n')

ITOL.write(annotation)
ITOL.close()
