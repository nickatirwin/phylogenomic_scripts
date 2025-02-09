###  PARSE HMMSCAN OUTPUT
###
###  10-03-2021
#
# Useful for making ITOL domain annotation files.
#
# usage: python parse_HMMscan_alignments.py proteinA.hmmscan 

import sys


# load in files and output
hmmfile = open(sys.argv[1],'r').readlines()

out = open(sys.argv[1]+'.parsed','w')
out.write('hmm\tdomain\tseq\tscore\tbias\tc-Evalue\ti-Evalue\thmmfrom\thmm_to\talifrom\tali_to\tenvfrom\tenvto\tacc\n')

hmm = hmmfile[6].split(' ')[-1].strip()

# parse hmmscan file
seq_d = {}

n = 0
while n < len(hmmfile):
    if hmmfile[n].startswith('Query:'):
        seq = hmmfile[n].split(' ')[-3].strip()
        n += 1
    try:
        if seq not in seq_d.keys():
            seq_d[seq] = []
    except:
        pass
    if hmmfile[n].startswith('>>'):
        domain = hmmfile[n].split(' ')[1].strip()
        n += 1
    elif ('!' in hmmfile[n]) or ('?' in hmmfile[n]):
        m = 1    
        while ('!' in hmmfile[n]) or ('?' in hmmfile[n]):
            if '!' in hmmfile[n]:
                info = hmmfile[n].split('!')[1]
            elif '?' in hmmfile[n]:
                info = hmmfile[n].split('?')[1]
            while "  " in info:
                info = info.replace("  ", " ")
            info = info.replace(" ", '\t').strip().split('\t')
            seq_d[seq].append([hmm,domain,seq+'_'+str(m),info[0],info[1],info[2],info[3],info[4],info[5],info[7],info[8],info[10],info[11],info[13]])
            n += 1
            m += 1
    else:
        n += 1

# filter out overlapping domains (if >50% overlapping, take the one with the better conditional e-value)

# set the overlap threshold - set to 0 if you want to keep all domains regardless of overlap
overlap = 0.5

for seq in seq_d:
    rm = []
    for domain in seq_d[seq]:
        for domain2 in seq_d[seq]:
            if domain != domain2:
                if (int(domain2[11]) >= int(domain[11])):
                    if len(set(range(int(domain2[11]),int(domain2[12])))-set(range(int(domain[11]),int(domain[12]))))/len(set(range(int(domain2[11]),int(domain2[12])))) <= overlap:
                        if (float(domain2[5]) > float(domain[5])):
                            rm.append(domain2)
                        else:
                            rm.append(domain)
    n = 0
    keepers = []
    while n < len(seq_d[seq]):
        if seq_d[seq][n] not in rm:
            keepers.append(seq_d[seq][n])
        n += 1
    seq_d[seq] = keepers

# output the results         
sep = '\t'       
for seq in seq_d:
    for domain in seq_d[seq]:
        out.write(sep.join(domain)+'\n')
out.close()
        



