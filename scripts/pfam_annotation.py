# ANNOTATE A TREE WITH PFAM DOMAINS
#
# Annotate sequences within a tree using Pfam annotations
#
# Requires a Pfam hmmscan file
#

'''
# download the Pfam database
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

# unzip the database
gunzip Pfam-A.hmm.gz

# compress the database for use with the HMMscan function
hmmpress Pfam-A.hmm

# annotate your proteins
hmmscan -E 1e-5 --incE 1e-5 --domE 1e-5 --cpu 5 -o [fasta file].hmmscan Pfam-A.hmm [fasta file]

# run the script
python pfam_annotation.py proteinA.fasta proteinA.fasta.hmmscan
'''

# usage: python pfam_annotation.py [fasta file] [hmmscan file]
#
# Outputs an annotation file ('.pfam.annotation') and a parsed hmmscan file ('.hmmscan.parsed')

import sys

# first parse the hmmscan file

# load in files and output
hmmfile = open(sys.argv[2],'r').readlines()

out = open(sys.argv[1]+'.hmmscan.parsed','w')
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
    elif '!' in hmmfile[n]:
        m = 1
        while '!' in hmmfile[n]:
            info = hmmfile[n].split('!')[1]
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
                if (int(domain2[9]) >= int(domain[9])):
                    if len(set(range(int(domain2[9]),int(domain2[10])))-set(range(int(domain[9]),int(domain[10]))))/len(set(range(int(domain2[9]),int(domain2[10])))) <= overlap:
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

# make the figtree annotation file

fasta = open(sys.argv[1],'r').readlines()

fasta_seqs = []
for s in fasta:
    if s.startswith('>'):
        fasta_seqs.append(s.split('>')[1].strip())

pfam = open(sys.argv[1]+'.hmmscan.parsed','r').readlines()

pfam_d = {}
current_seq = pfam[1].split('\t')[2].rsplit('_',1)[0]
s_d = {}
sep = '+'
for s in pfam[1:]:
    if s.split('\t')[2].rsplit('_',1)[0] == current_seq:
        s_d[int(s.split('\t')[9])] = s.split('\t')[1]
    else:
        starts = list(s_d.keys())
        starts.sort()
        domains = []
        for n in starts:
            domains.append(s_d[n])
        pfam_d[current_seq] = sep.join(domains)
        s_d = {}
        s_d[int(s.split('\t')[9])] = s.split('\t')[1]
        current_seq = s.split('\t')[2].rsplit('_',1)[0]

out = open(sys.argv[1]+'.pfam.annotation','w')
out.write('seq\tprotein_domain\n')
for s in fasta_seqs:
    try:
        out.write(s+'\t'+pfam_d[s]+'\n')
    except:
        out.write(s+'\tNoDomains\n')
out.close()
        
        

