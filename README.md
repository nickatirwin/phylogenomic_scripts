# Phylogenomic scripts

A series of scripts that are useful for analyzing phylogenomic data.
The descriptions of the scripts have been split into a few different sections.

Note: to get the most out of these scripts, it is worth formatting the headers of all fasta files as taxaID.proteinID (e.g., 9606.Q6NXT2). NCBI taxonomy information can be found at https://www.ncbi.nlm.nih.gov/taxonomy.

## Installation
Install dependencies:
```
conda install --yes --file requirements.txt
# or
mamba install --yes --file requirements.txt
```

## Database assembly

### download_reference_proteomes.py
Download a taxonomically balanced set of proteomes from UniProt Reference Proteomes (https://www.uniprot.org/proteomes). From a Uniprot proteome list, the script will download the best proteome per rank1, up to a max of n proteomes per family rank 2, but will include mandatory taxa, can exclude certain taxa, and will assess quality using busco (or protein number) given minimum thresholds.
The input is a uniprot table in the format Entry,Organism,Organism ID, protein count, BUSCO, CPD, Genome Assembly ID.
![alt text](misc/uniprot_image.png)
Note that for cetain taxa (like viruses) protein count may be a better quality metric than busco.

```
Options:
  -h, --help          show this help message and exit
  -i --input          Table of UniProt proteomes
  -r1 --rank1         Lower specified rank (best proteomer per rank1...)
  -n --max_proteomes  Maximum proteomes per upper rank (...up to n proteomes per rank2)
  -r2 --rank2         Upper specified rank
  -e --exclude        List of taxa ids to exclude
  -m --mandatory      Mandatory taxa to keep
  -q --quality        Quality metric to use (busco or counts (number of sequences)) [default = busco]
  -b --busco          Minimum busco to keep [default = 0]
  -s --sequences      Minimum number of sequences to keep [default = 0]
  -a --all            Download all proteomes from tsv file [default = no]

# example usage:
python download_reference_proteomes.py -i uniprot-fungi.filt.tsv -r1 species -n 1 -r2 order -e 1813822,1603295,6029,56615 -m 559292,237561,765915,578462,1340429,1220926,44442,994334,763407,658196 -q busco -b 25 -s 250
```

### genbank_download.py
Given a list of genbank accessions, download the sequences, write them to a fasta file, and make a Figtree annotation file with taxonomy and metadata. The query can be included in the output file if desired.
```
python genbank_download.py [list of accessions] [query sequences] [nucleotide/protien]
# eg.
python genbank_download.py accession_list.txt query_sequences.fasta nucleotide
# eg. without the query
python genbank_download.py accession_list.txt NA nucleotide
```

### genome_download.py
Download NCBI proteomes or genomes for a given taxonomic group or list of taxa.
```
python genome_download.py [list of taxids] [format - eg. protein OR genome]
# eg. Download all ciliate genomes
python genome_download.py 5878 genome
# eg. Download proteomes from a list of taxa IDs (each line is a different taxa ID)
python genome_download.py taxa.list protein
```

## Sequence and alignment processing
Scripts for processing sequence and alignment files

### filter_sequence_lengths.py
Filter out sequences from a fasta file that are above or below a certain size.
```
python filter_sequence_lengths.py [fasta] [min/max] [min_seq_length]
# eg.
python filter_sequence_lengths.py proteinA.fasta max 100
python filter_sequence_lengths.py proteinA.fasta min 50 100
```

### filter_alignment_by_percent_positions.py
Filter an alignment to remove sequences with a minimum number of sites. E.g., remove sequences that are 50% gaps.
```
python filter_alignment_by_percent_positions.py [alignment.fasta] [proportion]
# eg.
python filter_alignment_by_percent_positions.py proteinA.fasta.aln 0.2
```

### parse_HMMscan_alignments.py
Parse the output of HMMscan which is useful for looking at domain annotations or extracting domains. For domains that are overlapping by >50%, we take the domain with the better conditional e-value. To take all domains, replace the overlap value on line 55.
```
python parse_HMMscan_alignments.py [hmmscan output]
```
Note: the hmmscan file should be produced in the default output using a command such as:
```
hmmscan -E 1e-5 --incE 1e-5 --domE 1e-5 --cpu 5 -o proteinA.hmmscan Pfam-A.hmm proteinA.fasta
```

### domain_split.py
Split a fasta file into individual domains annotated using HMMscan. The input should include a parsed HMMscan output (see above).  
```
python domain_split.py [fasta file] [parsed hmmscan output]
# eg.
python domain_split.py proteinA.fasta proteinA.fasta.hmmscan.parsed
```

### remove_duplicates.py
Remove sequences that have identical headers from a fasta file. The first sequence with a unique header will be kept.
```
python remove_duplicates.py [fasta]
```

### rename_duplicates.py
Rename sequences that have identical headers from a fasta file.
```
python remove_duplicates.py [fasta]
```

### remove_illegal.py
Filter out problematic amino acids ('U','O','J','Z','B') and replace them with 'X'.
```
python remove_illegal.py [fasta]
```

### single_line_fasta.py
Turn a multiline fasta into a single line fasta (ie. where each sequence has two lines only)
```
python single_line_fasta.py [fasta]
```

### taxid_to_names.py
Translate taxonomy ids into species names in a fasta file. E.g., 9606.ABC1 turns into Homo_sapiens.ABC1. Note that the header form must be taxid.proteinid.
```
python taxid_to_names.py [fasta]
```

## Phylogenomics
Scripts for generating phylogenomic matrices, parsing the data, and recoding alignments

### concatenation.py

Concatenate a set of a alignments into a supermatrix for phylogenomic analyses. Each alignment should have a maximum of one sequence per species and species names should be denoted at the start of the headers (seperated by a period - e.g., >Homo_sapiens.proteinID).

The output is a concatenated alignment (.fasta) and a statistics files noting the percentage of genes and sites present in each species (.species_stats.tab) and the percentage of species with each gene (.gene_stats.tab).

```
python concatenation.py [list of alignment files]
# eg.
python concatenation.py '*.fasta.aln'
```

### dayhoff_recoding.py

Recode an amino acid alignment using Dayhoff groups. This can be useful for dealing with or assessing saturation and compositional heterogeneity (e.g., see Susko & Roger 2007, MBE, https://doi.org/10.1093/molbev/msm144). The script can recode an alignment using 4-state or 6-state dayhoff groups.

```
python dayhoff_recoding.py [alignment file] [4 or 6]
```

### fast_site_removal.py

Remove a given proportion of the fastest evolving sites in an alignment. Site rates should be estimated using IQ-Tree (https://github.com/Cibiv/IQ-TREE) using the -wsr option (produces a .rate file). Removal of fast evolving sites can be useful for assessing long branch attraction and assessing the phylogenetic support for a given topology.

```
python fast_site_removal.py [alignment file] [.rate file] [proportion of sites to remove]
```

### monophyly_constraint.py
Create a monophyly constraint file for doing topology testing. Colour the taxa that should be monophyletic in FigTree, save the tree, generate a monophyly constraint file, and re-run the constrained topology.
```
python monophyly_constraint.py [treefile_with_monophyletic otus coloured] [colour]
# eg
python monophyly_constraint.py treefile.coloured 0000ff
```

### exonerate_protein_prediction.py

Predict proteins from a genome using a given protein or proteome for reference (e.g. proteome of a closely related species, transcriptome, or a set of genes). The gene models are usually rough but may provide a useful starting place. Useful for estimating genome completeness or identifying missing genes for phylogenomic analyses. 

In brief, proteins are mapped to the genome using tBLASTn (e < 1e-5), mapped regions are extracted, query proteins are used as a model for Exonerate to identify exons, exons are combined into coding regions, and coding regions are translated. The resulting protein predictions are compared against the original proteome and protein models with a hit back to the original dataset are retained and clustered at 99% using Cd-hit to reduce redundancy.

```
python exonerate_protein_prediction.py [proteome.fasta] [genome.fasta] [threads] [genetic code (integer)]
# eg.
python exonerate_protein_prediction.py reference_proteome.fasta genome_scaffolds.fasta 20 6
```

## ITOL annotation
Scripts for automatically generating ITOL annotation files.

### ITOL_binary.py
Generate a presence absence distribution based on one or more fasta files.
```
python ITOL_binary.py '*fasta'
```

### ITOL_domain_annotation.py
Plot domain structures across a phylogeny. Domain annotations are based on parsed HMMScan annotations (see parse_HMMscan_alignments.py).
```
python ITOL_domain_annotation.py [fasta file] [parsed hmmscan output]
# eg.
python ITOL_domain_annotation.py proteinA.fasta ProteinA.fasta.hmmscan.parsed
```

### ITOL_heatmap.py
Generate a presence absence heatmap across a phylogeny using one or more fasta files.
```
python ITOL_heatmap.py '*fasta'
```

## FigTree annotation
Scripts for annotating and working with trees using FigTree (http://tree.bio.ed.ac.uk/software/figtree/)

### remove_and_extract.py

Remove or extract sequences from a fasta file based on leaf colouring in FigTree.

```
# remove red sequences (ff0000)
python remove_and_extract.py proteinA.fasta proteinA.fasta.treefile.coloured ff0000

# extract blue sequences (0000ff)
python remove_and_extract.py proteinA.fasta proteinA.fasta.treefile.coloured NA 0000ff

# remove red sequences (ff0000) and extract blue sequences (0000ff)
python remove_and_extract.py proteinA.fasta proteinA.fasta.treefile.coloured ff0000 0000ff
```

### blast_annotation.py

Annotate the sequences in a tree visualized with FigTree using best blast hits from SwissProt. This is useful for adding functional annotations to a tree.

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

Annotate  sequences in a tree visualized in FigTree with Pfam domains (or any annotations made using HMMER hmmscan). This is useful for adding functional annotations to a tree. This script requires an hmmscan output file. The script will parse the hmmscan file such that if two domains are overlapping by over 50%, the domain with the better conditional e-value will be taken.

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

Annotate sequences with extra taxonomic information (domain, supergroup, and species name). It uses NCBI Taxonomy and requires sequence names to be in the format taxaID.proteinID (e.g. 9606.Q53XC5 - where 9606 is Homo sapiens and Q53XC5 is a uniprot protein accession). The script also requires ete3 (https://github.com/etetoolkit/ete) and NCBI Taxonomy database.

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


