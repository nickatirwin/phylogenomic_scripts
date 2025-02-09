# GENBANK DOWNLOAD AND ANNOTATE
#
# given a list of genbank accessions, download the sequences, write them to a fasta file,
# and make a Figtree annotation file with taxonomy and metadata
#
# usage: python genbank_download.py [list of accessions] [query sequences] [nucleotide/protien]
# 
# eg. python genbank_download.py accession_list.txt query_sequences.fasta nucleotide
# or to skip query sequences
# eg. python genbank_download.py accession_list.txt NA nucleotide

import os
import sys
from Bio import Entrez
from ete3 import NCBITaxa

ncbi = NCBITaxa()

def reverse_d(d):
    rd = {}
    for k in d.keys():
        rd[d[k]] = k
    return rd

def taxonomy_classifier(taxid):
    d ={'domain':'','kingdom':'','phylum':'','class':'','order':'','family':'','genus':'','species':''}
    # try to get the desired level but if not take family and then genus
    taxonomy = reverse_d(ncbi.get_rank(ncbi.get_lineage(int(taxid))))
    # record information
    lineage = ncbi.get_lineage(int(taxid))
    try:
        d['domain'] = lineage[2]
    except:
        d['domain'] = 'NA'
    # check for viruses
    if 10239 in lineage:
        d['domain'] = '10239'
    try:
        d['kingdom'] = taxonomy['kingdom']
    except:
        try:
            d['kingdom'] = taxonomy['subkingdom']
        except:
            try:
                d['kingdom'] = lineage[3]
            except:
                d['kingdom'] = 'NA'
    try:
        d['phylum'] = taxonomy['phylum']
    except:
        try:
            d['phylum'] = taxonomy['superphylum']
        except:
            try:
                d['phylum'] = taxonomy['subphylum']
            except:
                d['phylum'] = 'NA'
    try:
        d['class'] = taxonomy['class']
    except:
        try:
            d['class'] = taxonomy['superclass']
        except:
            try:
                d['class'] = taxonomy['subclass']
            except:
                d['class'] = 'NA'               
    try:
        d['order'] = taxonomy['order']
    except:
        try:
            d['order'] = taxonomy['superorder']
        except:
            try:
                d['order'] = taxonomy['suborder']
            except:
                try:
                    d['order'] = taxonomy['infraorder']
                except:
                    try:
                        d['order'] = taxonomy['parvorder']
                    except:
                        d['order'] = 'NA'
    try:
        d['family'] = taxonomy['family']
    except:
        try:
            d['family'] = taxonomy['superfamily']
        except:
            try:
                d['family'] = taxonomy['subfamily']
            except:
                d['family'] = 'NA'
    try:
        d['genus'] = taxonomy['genus']
    except:
        d['genus'] = 'NA'
    try:
        d['species'] = taxonomy['species']
    except:
        d['species'] = 'NA'
    # if all are 'NAs' then assign as environmental
    if len(list(set([d[x] for x in d.keys()]))) == 1:
        for t in d:
            d[t] = 'environmental'
    # make it readable
    for t in d:
        if (d[t] != 'NA') and (d[t] != 'environmental'):
            d[t] = ncbi.translate_to_names([d[t]])[0]    
    return(d)

# load accessions
accs = open(sys.argv[1],'r').readlines()
accs = list(set([a.strip() for a in accs]))

# download the sequences
Entrez.email = 'sample@gmail.com'
gb_file = sys.argv[1]+".gbk"
print("Downloading...")
if os.path.isfile(gb_file):
    og_gb_file = gb_file
    n = 1
    while os.path.isfile(gb_file):
        gb_file = og_gb_file+str(n)
        n += 1

seq_type = sys.argv[3]
seq_type = seq_type.lower()
if not os.path.isfile(gb_file):
    with Entrez.efetch(db=seq_type, id=accs, rettype="gb", retmode="text") as net_handle: 
        with open(gb_file, "w") as out_handle:
            out_handle.write(net_handle.read())
    print("Saved.")
    
# read and parse sequences and get metadata
seq_d = {} # accession: sequence, domain, supergroup, kingdom, phylum, class, order, genus, family, species, isolation_source, country, taxaid 
gb_file = open(sys.argv[1]+'.gbk','r').read().split('//\n')[:-1]
for seq in gb_file:
    # get the accession
    accession = seq.split('ACCESSION')[1].split('\n')[0].strip().replace(' ','_')
    # get the sequence
    sequence = seq.split('ORIGIN')[1].strip().replace(' ','').split('\n')
    sequence = str('').join([s.strip('1234567890') for s in sequence])
    sequence.replace('\n','')
    # get the taxonomy
    sp = seq.split('ORGANISM')[1].split('\n')[0].strip()
    try:
        taxid = ncbi.get_name_translator([sp])[sp][0]
    except:
        try:
            sp2 = sp.rsplit(' ',1)[0]
            taxid = ncbi.get_name_translator([sp2])[sp2][0]
        except:
            try:
                sp3 = sp.split(' ')[0]
                taxid = ncbi.get_name_translator([sp3])[sp3][0]
            except:
                taxid = 1
    t = taxonomy_classifier(taxid)
    t['taxaid'] = str(taxid)
    # get metadata
    metadata = seq.split('FEATURES')[1].split('ORIGIN')[0].split('/')
    metadata_d = {}
    for m in metadata:
        if '=' in m:
            m = m.strip().split('=')
            metadata_d[m[0]] = str(' ').join([x.strip() for x in m[1].rsplit('"',1)[0].replace('"','').split('\n')])
    if 'isolation_source' not in metadata_d.keys():
        metadata_d['isolation_source'] = 'NA'
    if 'country' not in metadata_d.keys():
        metadata_d['country'] = 'NA'
    # record data
    seq_d[accession] = [sequence,t['taxaid'],t['domain'],t['kingdom'],t['phylum'],t['class'],t['order'],t['family'],t['genus'],t['species'],metadata_d['isolation_source'],metadata_d['country']]

# add in query sequencesi
if sys.argv[2] != 'NA':
    fasta = open(sys.argv[2],'r').read().split('>')[1:]
    for seq in fasta:
        name = seq.split('\n')[0].split(' ')[0]
        seq_d[name] = [seq.split('\n',1)[1].replace('\n',''),seq.split('.')[0],'Query','Query','Query','Query','Query','Query','Query',name,'Query','Query']

# output
out = open(sys.argv[1]+'.fasta','w')
for seq in seq_d:
    out.write('>'+seq_d[seq][1]+'.'+seq+'\n'+seq_d[seq][0]+'\n')
out.close()

out = open(sys.argv[1]+'.annotation','w')
out.write('accession\ttaxaid\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tisolation_source\tcountry\n')
for seq in seq_d:
    out.write(seq_d[seq][1]+'.'+seq + '\t'+ str('\t').join(seq_d[seq][1:])+'\n')
out.close()
