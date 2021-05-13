# Figtree-Analysis

A series of scripts for generating annotation files and manipulating fasta files with FigTree.

## Overview

These annotation files are useful for examining species trees, identifying orthologs, and cleaning datasets for phylogenomics. FigTree is a graphical tree viewer available at https://github.com/rambaut/figtree/releases.

## Annotation scripts

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


