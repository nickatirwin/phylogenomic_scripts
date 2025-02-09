# download_reference_proteomes.py
#
# required: ncbi-datasets-cli (conda install -c conda-forge ncbi-datasets-cli)
#
# usage: python download_reference_proteomes.py [uniprot_download_list] [rank1] [max proteomes per rank2] [rank2] [taxa to exclude] [mandatory taxa]
#        [scoring mode (busco/protein count)] [minimum busco] [minimum proteins] 
#
# eg. from a uniprot proteome list, download the best proteome per genus (rank1), up to a max of n proteomes per family (rank 2), excluding these taxa,
#     including these taxa, and assessing proteome quality using busco with a minimum busco of 60% and 100 proteins 
#
# NOTE: the input uniprot table should be Entry,Organism,Organism ID, protein count, BUSCO, CPD, Genome Assembly ID

# load modules
import sys
import argparse
import subprocess
import random
from datetime import date
from glob import glob
from ete3 import NCBITaxa

ncbi = NCBITaxa()
print()

# load in arguments
# define arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", metavar='--input',help="Table of UniProt proteomes")
parser.add_argument("-r1", metavar='--rank1', help="Lower specified rank (best proteomer per rank1...)")
parser.add_argument("-n", metavar='--max_proteomes', help="Maximum proteomes per upper rank (...up to n proteomes per rank2)")
parser.add_argument("-r2", metavar='--rank2', help="Upper specified rank")
parser.add_argument("-e", metavar='--exclude', help="List of taxa ids to exclude")
parser.add_argument("-m", metavar='--mandatory',help="Mandatory taxa to keep")
parser.add_argument("-q", metavar='--quality', help="Quality metric to use (busco or counts (number of sequences)) [default = busco]")
parser.add_argument("-b", metavar='--busco',help="Minimum busco to keep [default = 0]")
parser.add_argument("-s", metavar='--sequences', help="Minimum number of sequences to keep [default = 0]")
parser.add_argument("-a", metavar='--all', help="Download all proteomes from tsv file [default = no]")
args = parser.parse_args()

# check for required inputs
if (args.i == None) or (args.r1 == None) or (args.n == None) or (args.r2 == None):
    print('\nError: Missing required inputs (-i, -r1, -n, -r2)\n')
    parser.print_usage()
    print()
    sys.exit()

# set default inputs
if args.e == None:
    args.e = 'NA'
if args.m == None:
    args.m = 'NA'
if args.q == None:
    args.q = 'busco'
if args.b == None:
    args.b = '0'
if args.s == None:
    args.s = '0'
if args.a == None:
    args.a = 'no'

proteome_list = args.i
rank1 = args.r1
maximum = int(args.n)
rank2 = args.r2

exclude = args.e
exclude = exclude.split(',')

mandatory = args.m
mandatory = mandatory.split(',')

mode = args.q
minimum_busco = float(args.b)
minimum_proteins = float(args.s)

# define functions
# reverse dictionary
def reverse_d(d):
    rd = {}
    for k in d.keys():
        rd[d[k]] = k
    return rd
# map back taxonomy
levels = ['superkingdom','kingdom','subkingdom','superphylum','phylum','subphylum','superclass','class','subclass','superorder','order','suborder','infraorder','parvorder','superfamily','family','subfamily','genus','species']
def map_back_taxonomy(taxid,level):
    # try to get the desired level but if not take the nearest rank
    try:
        out = reverse_d(ncbi.get_rank(ncbi.get_lineage(int(taxid))))[level]
    except:
        try:
            available_ranks = list(reverse_d(ncbi.get_rank(ncbi.get_lineage(int(taxid)))).keys())
            found_ranks = []
            for i in levels:
                if i in available_ranks:
                    found_ranks.append(i)
                else:
                    found_ranks.append('')
            i = levels.index(level)
            rank, found, n = '', False, 0
            while (found == False) and (n < len(levels)):
                if found_ranks[i-n] != '':
                    found = True
                    rank = i-n
                    break
                else:
                    if found_ranks[i+n] != '':
                        found = True
                        rank = i+n
                        break
                n += 1
            if rank != '':
                out = reverse_d(ncbi.get_rank(ncbi.get_lineage(int(taxid))))[found_ranks[rank]]
            else:
                out = reverse_d(ncbi.get_rank(ncbi.get_lineage(int(taxid))))['genus']
        except:
                out = 'NA'
                print('Error: taxonomy could not be mapped back for: '+str(taxid))
    return(str(out))
    
# clean sequences of illegal characters
def sequence_cleaner(s):
    clean_s = s.replace('*','').replace('U','X').replace('B','X').replace('O','X').replace('J','X').replace('Z','X')
    return clean_s
    
# load in the proteomes and record metadata
dataset_d = {} # proteome_ID: [taxaid, rank1, rank2, protein_count, busco, assembly_ID, mandatory] 
uniprot = open(proteome_list,'r').readlines()[1:]
print('Loading datasets...')
for p in uniprot:
    # get proteome ID
    pid = p.split('\t')[0]
    # get the taxaid
    taxaid = p.split('\t')[2]
    # if taxaid cannot be read, exclude
    try:
        lineage = ncbi.get_lineage(taxaid)
    except:
        print('Error. Could not read taxa ID: ' + taxaid)
        continue
    # check whether the taxa should be excluded
    skip = 'no'
    if exclude != ['NA']:
        for t in lineage:
            if str(t) in exclude:
                skip = 'yes'
    # check whether the taxa is mandatory
    m = 'no'
    if mandatory != ['NA']:
        for t in lineage:
            if str(t) in mandatory:
                m = 'yes'
    # record rank1 and rank2
    r1 = map_back_taxonomy(taxaid,rank1)
    r2 = map_back_taxonomy(taxaid,rank2)
    # record sequence count
    protein_count = float(p.split('\t')[3])
    # record percentage of present busco
    if ',M:' in p: # only check proteomes with BUSCO
        busco = 100-float(p.split(',M:')[1].split('%')[0].strip())
    else:
        busco = 0
    # get assembly ID
    assembly = p.split('\t')[6].strip()
    # record data
    if skip == 'no':
        if (busco >= minimum_busco) and (protein_count >= minimum_proteins):
            dataset_d[pid] = [taxaid,r1,r2,protein_count,busco,assembly,m]

# select the best proteome per rank1
best_r1 = {} # rank1_id: proteome 
for p in dataset_d:    
    try:
        if mode == 'busco':
            if dataset_d[best_r1[dataset_d[p][1]]][4] < dataset_d[p][4]:
                best_r1[dataset_d[p][1]] = p
        elif mode == 'counts':
            if dataset_d[best_r1[dataset_d[p][1]]][3] < dataset_d[p][3]:
                best_r1[dataset_d[p][1]] = p
    except:
        best_r1[dataset_d[p][1]] = p

# select the best proteomes per rank2 (to a maximum of n proteomes)
# first record all of the proteomes from the best of rank1 and their scores for each rank2
r2_d = {} # {rank2: {proteome:score}} ie nested dictionary
for t in best_r1:
    p = best_r1[t]
    try:
        if mode == 'busco':
            r2_d[dataset_d[p][2]][p] = dataset_d[p][4]
        elif mode == 'counts':
            r2_d[dataset_d[p][2]][p] = dataset_d[p][3]
    except:
        if mode == 'busco':
            r2_d[dataset_d[p][2]] = {p:dataset_d[p][4]}
        elif mode == 'counts':
            r2_d[dataset_d[p][2]] = {p:dataset_d[p][3]}

# get the best n proteomes per rank2
best_r2 = {}
for r2 in r2_d:
    # get the dictionary for each rank2 and sort by score
    x = r2_d[r2]
    x = sorted(x, key=x.get)
    x.reverse()
    # take the best rank1 proteomes up to a maximum of n per rank2
    best_r2[r2] = x[0:maximum]  

# add in mandatory proteomes
for p in dataset_d:
    if dataset_d[p][6] == 'yes':
        mandatory_r1 = dataset_d[p][1]
        mandatory_r2 = dataset_d[p][2]
        # check the best rank2s, record mandatory proteomes and expendable ones
        mandatory_list,non_mandatory_list = [p],[]
        for p in best_r2[mandatory_r2]:
            if dataset_d[p][6] == 'yes':
                mandatory_list.append(p)
            elif dataset_d[p][2] == mandatory_r1:
                pass
            else:
                non_mandatory_list.append(p)
        mandatory_list = list(set(mandatory_list))
        non_mandatory_list = list(set(non_mandatory_list))
        final_list = mandatory_list   
        # if there are not enough mandatory proteomes to fill the quota, select additional proteomes
        final_list.extend(non_mandatory_list[0:maximum - len(final_list)])
        # replace the old list of proteomes to keep with the new list including the mandatory proteome
        best_r2[mandatory_r2] = final_list

# download the best proteomes and output a datafile and tree for the selected proteomes
print('\nDownloading reference proteomes...')
subprocess.call('mkdir downloads',shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# get a list of proteomes to download
proteomes = []
for r2 in best_r2:
    for p in best_r2[r2]:
        proteomes.append(p)

# if -a (-all) specified, download all proteomes
if args.a == 'yes':
    for p in dataset_d:
        if args.q == 'busco':
            if (dataset_d[p][4] >= float(args.b)):
                proteomes.append(p)
        if args.q == 'counts':
            if (dataset_d[p][3] >= float(args.s)):
                proteomes.append(p)
    proteomes = list(set(proteomes))


# download proteomes
pop, pn, total = [], 1, len(proteomes)
for p in proteomes:
    ### progress
    sys.stdout.write('\r')
    sys.stdout.write(str(pn) + '/' + str(total))
    sys.stdout.flush()
    pn += 1
    ###
    gca = dataset_d[p][5]
    skip = False
    # try to download the GCA
    try:
        subprocess.call('rm -r ncbi_dataset/',shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        subprocess.call('rm md5sum.txt ncbi_dataset.zip README.md',shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        subprocess.call('datasets download genome accession --include protein ' + gca,shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        subprocess.call('unzip ncbi_dataset.zip',shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        fasta = open('ncbi_dataset/data/'+gca+'/protein.faa','r').read().split('>')[1:]
    # if that fails, try to download the GCF
    except:
        try:
            # remove old files
            subprocess.call('rm -r ncbi_dataset/',shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            subprocess.call('rm md5sum.txt ncbi_dataset.zip README.md',shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            # download GCF
            gcf = gca.replace('GCA','GCF').split('.')[0]
            subprocess.call('datasets download genome accession --include protein ' + gcf,shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            subprocess.call('unzip ncbi_dataset.zip',shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            gcf = 'GCF_'+str(subprocess.check_output("find ncbi_dataset/ -name 'protein.faa'",shell=True, stderr=subprocess.STDOUT)).split('GCF_')[1].split('/')[0]
            fasta = open('ncbi_dataset/data/'+gcf+'/protein.faa','r').read().split('>')[1:]
            gca = gcf
            dataset_d[p][5] = gcf
        # if that fails, download direct from UniProt
        except:
            try:
                taxid = dataset_d[p][0]
                if 2759 in ncbi.get_lineage(taxid):
                    taxa = 'Eukaryota'
                elif 10239 in ncbi.get_lineage(taxid):
                    taxa = 'Viruses'
                elif 2157 in ncbi.get_lineage(taxid):
                    taxa = 'Archaea'
                else:
                    taxa = 'Bacteria'
                subprocess.call('wget --no-parent --no-passive -r https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/'+taxa+'/'+p+'/'+p+'_'+taxid+'.fasta.gz',shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                f = str(subprocess.check_output("find ftp.uniprot.org/ -name '*fasta.gz'",shell=True, stderr=subprocess.STDOUT)).split("'")[1].split("\\")[0]
                subprocess.call('gunzip ' + f,shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                f = f.replace('.gz','')
                fasta = open(f,'r').read().replace('sp|','').replace('tr|','').replace('|','.').replace('->','')
                fasta = fasta.split('>')[1:]
                gca = p
                subprocess.call('rm -r ftp.uniprot.org',shell=True)
            except:
                print(' Cannot download: ' + dataset_d[p][0])
                subprocess.call('rm -r ncbi_dataset/',shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                subprocess.call('rm md5sum.txt ncbi_dataset.zip README.md',shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                pop.append(p)
                skip = True
    if skip != True:
        out = open('downloads/'+gca.split('.')[0]+'_'+dataset_d[p][0]+'.renamed.fasta','w')
        for seq in fasta:
            out.write('>'+dataset_d[p][0]+'.'+seq.split('.')[0]+'\n'+sequence_cleaner(seq.split('\n',1)[1].replace('\n',''))+'\n')
        out.close()
        subprocess.call('rm -r ncbi_dataset/',shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        subprocess.call('rm md5sum.txt ncbi_dataset.zip README.md',shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# remove proteomes that couldn't be downloaded
for p in proteomes:
    if p in pop:
        proteomes.remove(p)

# output datafile
out = open(proteome_list.rsplit('.',1)[0]+'.selected.'+date.today().strftime("%y_%m_%d")+'.tab','w')
out.write('uniprot_id\tspecies\ttaxa_id\trank1\trank1_id\trank2\trank2_id\tproteins\tbusco\tassembly_id\n')
for p in proteomes:
    species = ncbi.translate_to_names([int(dataset_d[p][0])])[0]
    r1 = ncbi.translate_to_names([int(dataset_d[p][1])])[0]
    r2 = ncbi.translate_to_names([int(dataset_d[p][2])])[0]
    out.write(p+'\t'+species+'\t'+dataset_d[p][0]+'\t'+r1+'\t'+dataset_d[p][1]+'\t'+r2+'\t'+dataset_d[p][2]+'\t'+str(dataset_d[p][3])+'\t'+str(dataset_d[p][4])+'\t'+dataset_d[p][5]+'\n')
out.close()

# print NCBI taxonomy tree of species in the dataset
taxa = [int(dataset_d[p][0]) for p in proteomes]
t = ncbi.get_topology(taxa, intermediate_nodes=True)
t.convert_to_ultrametric()
for node in t.traverse():
    if not node.is_leaf():
        node.name = "INT"+node.name
t.write(format=1, outfile="NCBITaxonomy_tree."+proteome_list+".nwk")

print()
