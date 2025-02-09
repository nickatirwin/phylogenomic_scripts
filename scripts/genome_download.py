# download genomic data

# requires:
# conda install conda-forge::ncbi-datasets-cli

# usage: python genome_download.py [list of taxids] [format - eg. protein OR genome]

import sys
import subprocess
import pandas as pd

# load taxids to download (either from file or directly provided)
try:
    taxids = open(sys.argv[1],'r').readlines()
except:
    taxids = sys.argv[1].split(',')

# determine file end
if sys.argv[2] == 'protein':
    fe = 'faa'
else:
    fe = 'fna'

for t in taxids:
    t = t.strip()
    try:
        try:
            subprocess.call('datasets download genome taxon ' + str(t) + ' --reference --include ' + sys.argv[2] +' --filename ncbi_dataset.zip',shell=True)
            subprocess.call('unzip ncbi_dataset.zip',shell=True)
            fnames = str(subprocess.check_output('find ncbi_dataset/ -name "*'+fe+'"',shell=True)).strip().split("'")[1].split("\\n")[0:-1]
        except:
            try:
                subprocess.call('rm -r md5sum.txt README.md ncbi_dataset.zip ncbi_dataset/',shell=True)
                subprocess.call('datasets download genome taxon ' + str(t) + ' --include ' + sys.argv[2] + '  --filename ncbi_dataset.zip',shell=True)
                subprocess.call('unzip ncbi_dataset.zip',shell=True)
                fnames = str(subprocess.check_output('find ncbi_dataset/ -name "*'+fe+'"',shell=True)).strip().split("'")[1].split("\\n")[0:-1]
            except:
                new_t = str(ncbi.get_lineage(t)[-2])
                subprocess.call('datasets download genome taxon ' + str(new_t) + ' --reference --include ' + sys.argv[2] +' --filename ncbi_dataset.zip',shell=True)
                subprocess.call('unzip ncbi_dataset.zip',shell=True)
                fnames = str(subprocess.check_output('find ncbi_dataset/ -name "*'+fe+'"',shell=True)).strip().split("'")[1].split("\\n")[0:-1]
        # get taxonomy id from downloaded info
        metadata = pd.read_json('ncbi_dataset/data/assembly_data_report.jsonl', lines=True)
        # reformat and clean up downloaded data
        for fname in fnames:
            # get taxonomy id for each downloaded dataset
            try:
                t = str(metadata[metadata['accession']==fname.split('/')[2]]['assemblyInfo'].iloc[0]['biosample']['description']['organism']['taxId'])
            except:
                t = input('Enter taxid for ' + fname + ': ')
            # read and output genomic data
            out = open(t+'.'+fname.split('/')[2]+'.'+sys.argv[2]+'.fasta','w')
            fasta = open(fname,'r').read().split('>')[1:]
            for seq in fasta:
                out.write('>'+t+'.'+seq.split('\n',1)[0].split(' ')[0]+'\n'+seq.split('\n',1)[1].replace('\n','')+'\n')
            out.close()
        # remove temp files
        subprocess.call('rm -r md5sum.txt README.md ncbi_dataset.zip ncbi_dataset/',shell=True)
    except:
        print('\nERROR: '+str(t)+'\n')

