# fasta taxid to names
#
# convert the form >9606.ASDASC to >Homo_sapiens.ASDASC
#
# usage: python taxid_to_names.py myfasta.fasta

import sys
from ete3 import NCBITaxa

ncbi = NCBITaxa()

fasta = open(sys.argv[1],'r').readlines()
out =  open(sys.argv[1]+'.renamed','w')

for line in fasta:
    if line.startswith('>'):
        name = ncbi.translate_to_names([int(line.split('.')[0].strip('>'))])[0].replace(' ','_').replace('.','')
        out.write('>'+name+'.'+line.split('.',1)[1].strip()+'\n')
    else:
        out.write(line)
out.close()
        
