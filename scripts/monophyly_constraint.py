# MAKE CONSTRAINT FILE
#
# Make a monophyletic constraint for coloured taxa for use in topology testing and constrained tree searches
#
# USAGE: python monophyly_constraint.py [treefile_with_monophyletic group coloured] [colour]
# USAGE: python monophyly_constraint.py treefile.coloured 0000ff

import sys

treefile = open(sys.argv[1], 'r').readlines()[4:] # skip the begining strings
colour = '#' + sys.argv[2].strip('#')

other_taxa = []
monophyletic_taxa = []

for line in treefile:
    if line.strip() != ';':
        if colour in line:
            monophyletic_taxa.append(line.split('[')[0].strip().strip("'"))
        else:
            other_taxa.append(line.split('[')[0].strip().strip("'"))        
    elif line.strip() == ';':
        break

out = open(sys.argv[1]+'.monophyly_constraint','w')

monoclade = '('
for t in monophyletic_taxa:
    monoclade = monoclade + t + ','
monoclade = monoclade.strip(',') + ')'

constrained_tree = '('+monoclade+','
for t in other_taxa:
    constrained_tree = constrained_tree + t + ','
constrained_tree = constrained_tree.strip(',')+');'

out.write(constrained_tree)
out.close()


