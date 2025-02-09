# Phylogenomic scripts

A series of scripts that are useful for analyzing phylogenomic data.

## Database assembly scripts


## Scripts

### exonerate_protein_prediction.py

Predict proteins from a genome using a given proteome for reference (e.g. proteome of a closely related species, transcriptome, or a set of genes). The gene models can be rough but may provide a useful starting place and is useful for estimating genome completeness or identifying genes for phylogenomic analyses. 

In brief, proteins are mapped to the genome using tBLASTn (e < 1e-5), mapped regions are extracted, query proteins are used as a model for Exonerate to identify exons, exons are combined into coding regions, and coding regions are translated. The resulting protein predictions are compared against the original proteome and protein models with a hit back to the original dataset are retained and clustered at 99% to reduce redundancy. The pipeline is fairly quick and tends to run in 10-30 minutes with 60 threads for the datasets I've tested.

Required dependencies:
```
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
```
Usage:
```
python exonerate_protein_prediction.py [proteome.fasta] [genome.fasta] [threads] [genetic code (integer)]
```
Example:
```
python exonerate_protein_prediction.py reference_proteome.fasta genome_scaffolds.fasta 20 6
```
### fast_site_removal.py

This script will remove a given proportion of the fastest evolving sites in an alignment. Site rates should be estimated using IQ-Tree (https://github.com/Cibiv/IQ-TREE) using the -wsr option (produces a .rate file). Removal of fast evolving sites can be useful for assessing long branch attraction and assessing the phylogenetic support for a given topology.

Usage:
```
python fast_site_removal.py [alignment file] [.rate file] [proportion of sites to remove]
```
Example:
```
python fast_site_removal.py fasta.aln fasta.aln.rate 0.25
```
### concatenation.py

Concatenate a set of a alignments into a supermatrix for phylogenomic analyses. Each alignment should have a maximum of one sequence per species and species names should be denoted at the start of the headers (seperated by a period - e.g., >Homo_sapiens.proteinID).

The output is a concatenated alignment (.fasta) and a statistics files noting the percentage of genes and sites present in each species (.species_stats.tab) and the percentage of species with each gene (.gene_stats.tab).

Usage:
```
python concatenation.py [list of alignment files]
```
Example (important to include the quotes):
```
python concatenation.py '*.fasta.aln'
```
### dayhoff_recoding.py

Recode an amino acid alignment using Dayhoff groups. This can be useful for dealing with or assessing saturation and compositional heterogeneity (e.g., see Susko & Roger 2007, MBE, https://doi.org/10.1093/molbev/msm144). The script can recode an alignment using 4-state or 6-state dayhoff groups.

Usage:
```
python dayhoff_recoding.py [alignment file] [4 or 6]
```
Example:
```
python dayhoff_recoding.py fasta.aln 4
```

A series of scripts for generating annotation files and manipulating fasta files with FigTree.

## Overview

These annotation files are useful for examining species trees, identifying orthologs, and cleaning datasets for phylogenomics. FigTree is a graphical tree viewer available at https://github.com/rambaut/figtree/releases.

## Scripts

### remove_and_extract.py

This script allows one to remove or extract sequences from a fasta file based on leaf colouring in FigTree. E.g., Colour in contamination red and remove those sequences from the original fasta file.

```
# remove red sequences (ff0000)
python remove_and_extract.py proteinA.fasta proteinA.fasta.treefile.coloured ff0000

# extract blue sequences (0000ff)
python remove_and_extract.py proteinA.fasta proteinA.fasta.treefile.coloured NA 0000ff

# remove red sequences (ff0000) and extract blue sequences (0000ff)
python remove_and_extract.py proteinA.fasta proteinA.fasta.treefile.coloured ff0000 0000ff
```

### blast_annotation.py

This script allows one to annotate the sequences in a tree using best blast hits from Uniprot. This is useful for adding functional annotations to a tree. This script requires a blast output file (in a tab-delimited form - outfmt 6) 

```
# download swiss-prot and unzip
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# replace spaces with underscore so that the full annotation is included in the blast output
sed -i 's/ /_/g' uniprot_sprot.fasta

# make a blast database
makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot.fasta.db -dbtype prot

# run the blast
blastp -query proteinA.fasta -db uniprot_sprot.fasta.db -outfmt 6 -evalue 1e-5 -out proteinA.fasta.blastp

python blast_annotation.py proteinA.fasta proteinA.fasta.blastp
```
### pfam_annotation.py

This script allows one to annotate the sequences in a tree with Pfam domains (or any annotations made using HMMER hmmscan). This is useful for adding functional annotations to a tree. This script requires an hmmscan output file. The script will parse the hmmscan file such that if two domains are overlapping by over 50%, the domain with the better conditional e-value will be taken.

```
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
```
### ncbi_taxid_to_annotation.py

This script allows one to annotate sequences with extra taxonomic information (domain, supergroup, and species name). It uses NCBI Taxonomy and requires sequence names to be in the format taxaID.proteinID (e.g. 9606.Q53XC5 - where 9606 is Homo sapiens and Q53XC5 is a uniprot protein accession). The script also requires ete3 (https://github.com/etetoolkit/ete) and NCBI Taxonomy database.

```
# install ete3
conda install -c bioconda ete3

# download and update the NCBI Taxonomy database (run in python)
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

# run the script
python ncbi_taxid_to_annotation.py ProteinA.fasta
```

