# hmmsplit.py 
#
# Split a fasta file by the domains identified using hmmscan.
#
# Usage: python hmmsplit.py proteinA.fasta proteinA.fasta.hmmscan.parsed
#

import sys

fasta = open(sys.argv[1],'r').read().split('>')[1:]
hmmscan = open(sys.argv[2],'r').readlines()[1:]

# make sequence dictionary

seq_d = {}
for line in fasta:
    seq_d[line.split('\n')[0].strip('>')] = line.split('\n',1)[1].replace('\n','').strip()
 
# make an annotation dictionary

annot_d = {}
for line in hmmscan:
    m = 1
    seq = line.split('\t')[2].strip()
    if seq not in annot_d:
        annot_d[seq] = {}
    domain = line.split('\t')[1].strip()
    if domain in annot_d[seq]:
        domain = domain + '.' + str(m)
        m += 1
        annot_d[seq][domain] = [line.split('\t')[11],line.split('\t')[12],line.split('\t')[5]]
    else:
        annot_d[seq][domain] = [line.split('\t')[11],line.split('\t')[12],line.split('\t')[5]]
        
# filter annotations that are overlapping (if >75% overlapping, take the one with the better conditional e-value)
for seq in annot_d:
    rm = []
    for domain in annot_d[seq]:
        for domain2 in annot_d[seq]:
            if domain != domain2:
                if (int(annot_d[seq][domain2][0]) >= int(annot_d[seq][domain][0])):
                    if len(set(range(int(annot_d[seq][domain2][0]),int(annot_d[seq][domain2][1])))-set(range(int(annot_d[seq][domain][0]),int(annot_d[seq][domain][1]))))/len(set(range(int(annot_d[seq][domain2][0]),int(annot_d[seq][domain2][1])))) <= 0.5:
                        if (float(annot_d[seq][domain2][2]) > float(annot_d[seq][domain][2])):
                            rm.append(domain2)
                        else:
                            rm.append(domain)
    for d in list(set(rm)):
        annot_d[seq].pop(d)

# output 
for seq in annot_d:
    for domain in annot_d[seq]:
        out = open(domain.split('.')[0]+'.'+sys.argv[1],'a')
        out.write('>'+seq+'.'+domain.replace('.','_')+'\n')
        out.write(seq_d[seq.rsplit('_',1)[0]][int(annot_d[seq][domain][0])-1:int(annot_d[seq][domain][1])]+'\n')
        out.close()

