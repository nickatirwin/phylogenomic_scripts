# exonerate_protein_prediction.py 

'''
Predict proteins from a genome using exonerate and a given proteome as a query (eg. proteome from
a closely related genome(s) or a transcriptome(s)).

Protein predictions, coding regions, and a gff file are in results/.

usage: python exonerate_protein_prediction.py [proteins.fasta] [genome.fasta] [threads] [genetic code (int)]

dependencies:

# exonerate
conda install -c bioconda exonerate
# blast
conda install -c bioconda blast
# diamond
conda install -c bioconda diamond
# parallel
conda install -c conda-forge parallel
# cd-hit
conda install -c bioconda cd-hit
'''

import subprocess
from glob import glob
import sys
from statistics import median
from datetime import datetime

print('\nPredicting proteins for genome: ' + sys.argv[2] + '\nUsing proteins from: ' + sys.argv[1])
print('Start time: ' + str(datetime.now()))
subprocess.call('mkdir temp/', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# rename headers in protein and genome files to make them compatible with the script

def header_cleaner(s):
    clean_s = s.replace('*','').replace('|','_').replace(';','_').replace(',','_').replace('(','').replace(')','').replace('[','').replace(']','').replace('/','')
    return clean_s

proteins = open(sys.argv[1],'r').readlines()
out = open('temp/proteins.renamed.fasta','w')
for p in proteins:
    if p.startswith('>'):
        out.write(header_cleaner(p.split('\n')[0].split(' ')[0])+'\n')
    else:
        out.write(p.strip().replace('*','')+'\n')
out.close()

genome = open(sys.argv[2],'r').read().split('>')
out = open('temp/genome.renamed.fasta','w')
for s in genome[1:]:
    out.write('>'+header_cleaner(s.split('\n')[0].split(' ')[0])+'\n'+s.split('\n',1)[1].strip().replace('*','')+'\n')
out.close()

# map proteins to genome using tblastn

print('\nMapping proteins to genome using tBLASTn...')
# Make blast databases
subprocess.call('makeblastdb -in temp/genome.renamed.fasta -out temp/genome.renamed.fasta.db -dbtype nucl', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
# find homologs using tblastn - take the best hit (-max_target_seqs 1, -evalue 1e-5)
subprocess.call('tblastn -query temp/proteins.renamed.fasta -db temp/genome.renamed.fasta.db -max_target_seqs 1 -num_threads ' + sys.argv[3] + ' -evalue 1e-3 -outfmt 6 -out temp/genome_mapping.tblastn -db_gencode ' + sys.argv[4], shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# go through blast outputs and identify homology regions (take the max upstream and minumum downstream sites based on all HSPs)
print('Parsing blast output...')

# define the maximum gap size - eg. if multiple hits on one scaffold
max_gap_size = 2000

# load in blast files
blast = open('temp/genome_mapping.tblastn','r').readlines()
blast_d = {} # query_m: [scaffold, [start/stops], strand]

def strand_check(blast_line):
    if int(blast_line.split('\t')[8]) > int(blast_line.split('\t')[9]):
        strand = '-'
    else:
        strand = '+'
    return strand

def start_stop(blast_line, strand):
    if strand == '+':
        start = int(blast_line.split('\t')[8])
        stop = int(blast_line.split('\t')[9])
    else:
        start = int(blast_line.split('\t')[9])
        stop = int(blast_line.split('\t')[8])
    return [start, stop]

# parse blast file to identify homology regions and strand info
previous_start, previous_stop = 0, 0
previous_query = 'NA'



for line in blast:
    if (line.split('\t')[0] == previous_query):
        strand = strand_check(line)
        start, stop = start_stop(line,strand)
        if (abs(stop-previous_start) < max_gap_size) or (abs(start-previous_stop) < max_gap_size):
            try:
                blast_d[line.split('\t')[0]+'_1'][1] += [int(line.split('\t')[8]),int(line.split('\t')[9])]
            except:
                blast_d[line.split('\t')[0]+'_1'] = [line.split('\t')[1],[int(line.split('\t')[8]),int(line.split('\t')[9])],0,strand_check(line)]
            if start < previous_start:
                previous_start = start
            if stop > previous_stop:
                previous_stop = stop
    else:
        try:
            blast_d[line.split('\t')[0]+'_1'][1] += [int(line.split('\t')[8]),int(line.split('\t')[9])]
        except:
            blast_d[line.split('\t')[0]+'_1'] = [line.split('\t')[1],[int(line.split('\t')[8]),int(line.split('\t')[9])],0,strand_check(line)]
        previous_strand = strand_check(line)
        previous_start, previous_stop = start_stop(line,previous_strand)[0], start_stop(line,previous_strand)[1]
        previous_query = line.split('\t')[0]        
        
for s in blast_d:
    stop = max(blast_d[s][1])
    start = min(blast_d[s][1])
    blast_d[s][1] = start
    blast_d[s][2] = stop

# write out all homology regions and coresponding proteins for use with exonerate

# create sequences directory
subprocess.call('mkdir temp/sequences/', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# write fasta files for each protein sequence
print('Extracting sequences...')
proteins = open('temp/proteins.renamed.fasta','r').read().split('>')[1:]
protein_d = {}
for p in proteins:
    protein_d[p.split('\n')[0]] = p.split('\n',1)[1].replace('\n','')
for p in blast_d:
    out = open('temp/sequences/'+p+'.fasta','w')
    out.write('>'+p+'\n'+protein_d[p.rsplit('_',1)[0]]+'\n')
    out.close()

# write fasta files for each homology region (reverse complement if mapping to the negative strand)

def reverse_complement(seq):
    complement = ''
    for b in seq.upper()[::-1]:
        if b == 'A':
            complement += 'T'
        elif b == 'T':
            complement += 'A'
        elif b == 'C':
            complement += 'G'
        elif b == 'G':
            complement += 'C'
    return complement

genome = open('temp/genome.renamed.fasta','r').read().split('>')[1:]
genome_d = {}
homology_d = {}
for s in genome:
    genome_d[s.split('\n')[0]] = s.split('\n',1)[1].replace('\n','')
for i in blast_d:
    out = open('temp/sequences/'+i+'.homology.fasta','w')
    if blast_d[i][-1] == '+':
        out.write('>'+i+'_'+blast_d[i][0]+'\n'+genome_d[blast_d[i][0]][blast_d[i][1]-1:blast_d[i][2]]+'\n')
        homology_d[i] = genome_d[blast_d[i][0]][blast_d[i][1]-1:blast_d[i][2]]
    else: # if on the negative strand output the reverse complement
        out.write('>'+i+'_'+blast_d[i][0]+'\n'+reverse_complement(genome_d[blast_d[i][0]][blast_d[i][1]-1:blast_d[i][2]])+'\n')
        homology_d[i] = reverse_complement(genome_d[blast_d[i][0]][blast_d[i][1]-1:blast_d[i][2]])
    out.close()

# run exonerate - runs in parallel (1 job per thread) - score threshold = 50, intron parameters are generic and probably don't need to be changed (min size is crucial)

print('Using Exonerate to identify exons...')
subprocess.call("find temp/sequences/ -type f -name '*.homology.fasta' | awk -F '.homology.fasta' '{print $1}' | cut -f 3 -d '/' | parallel -j " + sys.argv[3] + " 'exonerate -m p2g --geneticcode " + sys.argv[4] + " --maxintron 20000 --minintron 0 --score 25 --bestn 1 --showtargetgff -q temp/sequences/{}.fasta -t temp/sequences/{}.homology.fasta | egrep -w exon > temp/sequences/{}.exon.gff'", shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# extract coding regions based on exonerate exons and write to CDS fasta in results/ - if exons overlap, join them

print('Extracting coding regions...')
subprocess.call('mkdir results/', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
out = open('results/'+sys.argv[2]+'.cds.fasta','w')

gffout = open('results/'+sys.argv[2]+'.gff','w')
gffout.write('## gff-version 3.2.1\n## assembly: ' + sys.argv[2] + '\n## template: ' + sys.argv[1].replace('.renamed','')+'\n')

n = 1
for seq in blast_d:
    gff = open('temp/sequences/'+seq+'.exon.gff','r').readlines()
    if len(gff) == 0:
        continue
    else:
        start = int(gff[0].split('\t')[3])
        stop =  int(gff[0].split('\t')[4])
        sequence = ''
        extract_d = {start:stop}
        extract_starts = [start]
        for hit in gff[1:]:
            if (int(hit.split('\t')[3]) < stop) and (int(hit.split('\t')[4]) > stop):
                extract_end = int(hit.split('\t')[4])
                extract_start = stop + 1
                stop = extract_end
                start = start
            elif (int(hit.split('\t')[3]) < start) and (int(hit.split('\t')[4]) > start):
                extract_start = int(hit.split('\t')[3])
                extract_end = start-1
                stop = stop
                start = extract_start
            else:
                extract_start = (int(hit.split('\t')[3]))
                extract_end = (int(hit.split('\t')[4]))
                if extract_start < start:
                    start = extract_start
                if extract_end > stop:
                    stop = extract_end
            extract_d[extract_start] = extract_end
            extract_starts.append(extract_start)
        extract_starts.sort()
        c = 1
        gene_start = 1
        gene_end = 1
        for s in extract_starts:
            sequence += homology_d[seq][s-1:extract_d[s]]
            # output gff annotation
            if blast_d[seq][3] == '+':
                gff_start = str(blast_d[seq][1]+s)
                gff_end = str(blast_d[seq][1]+extract_d[s])
                if str(c) == '1':
                    gffout.write(blast_d[seq][0]+'\t.\texon\t'+gff_start+'\t'+gff_end+'\t.\t'+blast_d[seq][3]+'\t.\tID=seq'+str(n)+'_exon'+str(c)+';Parent=seq'+str(n)+'\n')
                    gene_start = gff_start
                    gene_end = gff_end
                else:
                    gffout.write(blast_d[seq][0]+'\t.\texon\t'+gff_start+'\t'+gff_end+'\t.\t'+blast_d[seq][3]+'\t.\tID=seq'+str(n)+'_exon'+str(c)+';Parent=seq'+str(n)+'\n')
                    gene_end = gff_end
            else:
                gff_start = str(blast_d[seq][2]-s)
                gff_end = str(blast_d[seq][2]-extract_d[s])
                if str(c) == '1':
                    gffout.write(blast_d[seq][0]+'\t.\texon\t'+gff_end+'\t'+gff_start+'\t.\t'+blast_d[seq][3]+'\t.\tID=seq'+str(n)+'_exon'+str(len(extract_starts)-c+1)+';Parent=seq'+str(n)+'\n')
                    gene_start = gff_start
                    gene_end = gff_end
                else:
                    gffout.write(blast_d[seq][0]+'\t.\texon\t'+gff_end+'\t'+gff_start+'\t.\t'+blast_d[seq][3]+'\t.\tID=seq'+str(n)+'_exon'+str(len(extract_starts)-c+1)+';Parent=seq'+str(n)+'\n')
                    gene_end = gff_end
            c += 1
        if blast_d[seq][3] == '+':
            gffout.write(blast_d[seq][0]+'\t.\tgene\t'+str(gene_start)+'\t'+str(gene_end)+'\t.\t'+blast_d[seq][3]+'\t.\tID=seq'+str(n)+'\n')
        else:
            gffout.write(blast_d[seq][0]+'\t.\tgene\t'+str(gene_end)+'\t'+str(gene_start)+'\t.\t'+blast_d[seq][3]+'\t.\tID=seq'+str(n)+'\n')
        out.write('>seq'+str(n)+'.'+seq+'_'+blast_d[seq][0]+'\n'+sequence+'\n')
        # ouput to gff fileseq
        n += 1
out.close()
gffout.close()

# translate proteins using fastatranslate

print('Translating coding regions...')
subprocess.call('fastatranslate -f results/'+sys.argv[2]+'.cds.fasta --geneticcode ' + sys.argv[4] + ' > temp/translation.fasta',shell = True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# clean the gnarly headers

translations = open('temp/translation.fasta','r').read()
out = open('temp/translation.fasta','w')
out.write(translations.replace(' ','_').replace(':','_').replace('[','').replace(']','').replace('(','').replace(')','')) 
out.close()

# blast proteins to initial dataset to compare - we will only take proteins that had homologs in the query dataset - seems to reduce noise without compromising data

print('Comparing predicted proteins to original proteome...')
subprocess.call('diamond makedb --in temp/proteins.renamed.fasta --db temp/proteins.db', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
subprocess.call('diamond blastp -c1 --query temp/translation.fasta --max-target-seqs 1 --db temp/proteins.db --outfmt 6 --sensitive --evalue 1e-5 --out temp/translation.blastp --threads ' + sys.argv[3], shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# record blast results
blast = open('temp/translation.blastp', 'r').readlines()
protein_blast_d = {}
for line in blast:
    protein_blast_d[line.split('\t')[0]] = [line.split('\t')[1], float(line.split('\t')[-2])]

# record translated proteins
translation = open('temp/translation.fasta','r').read().split('>')[1:]
translation_d = {}
for seq in translation:
    translation_d[seq.split('\n')[0]] = seq.split('\n',1)[1].replace('\n','')

# identify the best translation per protein (longest and with best (lowest e-value) blast hit) - translations are done in all 6 frames
seq_d = {}
for p in translation_d:
    try:
        protein_blast_d[p]
        if (seq_d[p.split('.')[0]][1] == 'NA'):
            seq_d[p.split('.')[0]] = [translation_d[p],protein_blast_d[p][0],protein_blast_d[p][1],len(translation_d[p].replace('*',''))]
        elif (seq_d[p.split('.')[0]][2] > protein_blast_d[p][1]):
            seq_d[p.split('.')[0]] = [translation_d[p],protein_blast_d[p],protein_blast_d[p][1],len(translation_d[p].replace('*',''))]
    except:
        try:
            if (seq_d[p.split('.')[0]][1] != 'NA'):
                continue
            elif (seq_d[p.split('.')[0]][3] < translation[p][1]):
                seq_d[p.split('.')[0]] = [translation_d[p],protein_blast_d[p][0],protein_blast_d[p][1],len(translation_d[p].replace('*',''))]
        except:
            seq_d[p.split('.')[0]] = [translation_d[p],'NA','NA',len(translation_d[p].replace('*',''))]
    
# output the predicted proteins and cluster at 99% - reduces redundancy but this can be changed to 100% on line 274 (change '-c 0.99' to '-c 1.0')

print('Writing finished proteins...')
out = open('results/'+sys.argv[2]+'.proteins.fasta','w')
for seq in seq_d:
    sequence = max(seq_d[seq][0].split('*')) # split at stop codons - take the longest uninterupted sequence
    if seq_d[seq][1] != 'NA':
        out.write('>'+seq+'_length_'+str(len(sequence))+'\n'+sequence+'\n')
out.close()

# cluster proteins
print('Clustering resulting proteins at 95%...')
subprocess.call('cd-hit -i results/'+sys.argv[2]+'.proteins.fasta -c 0.95 -G 0 -aL 0.1 -o temp/cluster -T ' + sys.argv[3], shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
subprocess.call('mv temp/cluster results/'+sys.argv[2]+'.proteins.fasta',shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# finish
#subprocess.call('rm -r temp/sequences',shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
print('\nPredicted proteins and CDS in results directory\nFinish time: ' + str(datetime.now())+'\n')
