# REMOVE AND EXTRACT COLOURED TAXA FROM TREE
# 
# Given a tree, with undesired taxa coloured, this script will remove 
# those taxa from the original fasta file (outputs *.fasta.cleaned).
# This script is compatible with phylogenies made using RAxML, FastTree,
# and IQ-Tree.
# 
# You can also extract sequences with a certain colour using the fourth argument.
#
# Usage: python remove_and_extract.py [fasta_file] [coloured_tree] [colour_code_to_remove] [colour_code_to_extract]
# 
# Example 1: Remove red taxa and extract blue taxa
# python remove_and_extract.py example.fasta coloured.tre ff0000 0000ff
# 
# Example 2: Just remove red taxa
# python remove_and_extract.py example.fasta coloured.tre ff0000
#
# Example 3: Just extract blue taxa
# python remove_and_extract.py example.fasta coloured.tre NA ff0000

import sys

colour = sys.argv[3]

# do you want to remove sequences?
if colour != 'NA':
    output = open(sys.argv[1] + '.cleaned', 'w')

# do you want to extract sequences?
extract = 'no'
if len(sys.argv) > 4:
    extract = 'yes'
    extract_colour = sys.argv[4]
    ext_out = open(sys.argv[1]+'.extract','w')

# load in the fasta file
fasta_file = open(sys.argv[1], 'r').read().split('>')
fasta = []
for seq in fasta_file[1:]:
    fasta.append('>'+seq.split('\n')[0])
    sequence = ''
    for part in seq.split('\n')[1:]:
        sequence = sequence + part.strip()
    fasta.append(sequence)

# load in coloured tree file
tree = open(sys.argv[2], 'r').readlines()

# get lines with the coloured taxa info
coloured_taxa = []
extract_taxa = []

for lines in tree:
    if colour != 'NA':
        if colour in lines:
            coloured_taxa.append(lines.split("[")[0].strip("\t").strip("'").replace('@','_'))
    if extract == 'yes':
        if extract_colour in lines:
            extract_taxa.append(lines.split("[")[0].strip("\t").strip("'").replace('@','_'))

# remove and extract coloured taxa from the fasta file:
if colour != 'NA':
    n = 0
    while n < len(fasta):
        if fasta[n].startswith('>') and fasta[n].split(' ')[0].strip('>').strip('\n').replace('@','_') not in coloured_taxa:
            output.write(fasta[n].strip('\n') + '\n' + fasta[n+1].strip('\n') + '\n')
            n += 2
        elif fasta[n].startswith('>') and fasta[n].split(' ')[0].strip('>').strip('\n').replace('@','_') in extract_taxa:
            ext_out.write(fasta[n].strip('\n') + '\n' + fasta[n+1].strip('\n') + '\n')
            n += 2
        else:
            n += 2

if len(sys.argv) > 4:
    n = 0
    while n < len(fasta):
        if fasta[n].startswith('>') and fasta[n].split(' ')[0].strip('>').strip('\n').replace('@','_') in extract_taxa:
            ext_out.write(fasta[n].strip('\n') + '\n' + fasta[n+1].strip('\n') + '\n')
            n += 2
        else:
            n += 2

if colour != 'NA':
    output.close()
if len(sys.argv) > 4:
    ext_out.close()

