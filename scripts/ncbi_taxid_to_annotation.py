# ANNOTATE A PHYLOGENY WITH TAXONOMY INFORMATION
#
# Create an annotation file for FigTree from taxa ids in the headers of a fasta
# 
# NOTE: fasta file headers must be in the form "NCBITaxaID.ProteinID" (e.g. 9606.Q53XC5)
#
'''
# install ete3
conda install -c bioconda ete3

# download and update the NCBI Taxonomy database (run in python)
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

# run the script
python ncbi_taxid_to_annotation.py ProteinA.fasta
'''

from ete3 import NCBITaxa
ncbi = NCBITaxa()
import sys
import subprocess

#load in the fasta file
lines = open(sys.argv[1], 'r').readlines()
output = open(sys.argv[1] + '.taxonomy.annotation', 'w')
otu_d = {}

#extract the headers (OTUs) and taxa ids
for i in lines:
    if i.startswith('>'):
        otu_d[i.strip('>').strip()] = [i.split('.')[0].strip('>').strip()]

#get the names corresponding to each of the taxa ids using grep
for i in otu_d:
    taxa = i.split('.')[0]
    try: # add taxon name
        otu_d[i].append(ncbi.get_taxid_translator([int(taxa)])[int(taxa)])
    except:
        otu_d[i].append('NA')
    try: # add domain
        otu_d[i].append(str(ncbi.get_taxid_translator([ncbi.get_lineage(int(taxa))[3]])[ncbi.get_lineage(int(taxa))[3]]))
    except:
        otu_d[i].append('NA')
    try:
        domain = 'NA'
        if 2759 in ncbi.get_lineage(int(taxa)):
            domain = 'Eukaryota'
        elif 10239 in ncbi.get_lineage(int(taxa)):
            domain = 'Virus'
        elif 2 in ncbi.get_lineage(int(taxa)):
            domain = 'Bacteria'
        elif 2157 in ncbi.get_lineage(int(taxa)):
            domain = 'Archaea'
        else:
            domain = 'Virus'
        otu_d[i].append(domain)
    except:
        otu_d[i].append('NA')

#output the results
output.write('otu\ttaxaid\tspecies\tsupergroup\tdomain\n')
sep = '\t'
for otu in otu_d:
    output.write(otu+'\t'+sep.join(otu_d[otu])+'\n')
output.close()

