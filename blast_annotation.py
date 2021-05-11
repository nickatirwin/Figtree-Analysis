# ANNOTATE A TREE WITH BLAST HITS
#
# Annotate sequences within a tree using protein annotations from SWISS-PROT
#
# Requires a BLAST output file using SWISS-PROT and -outfmt 6
#
# To make this:
# 1. Download swiss-prot and unzip
#    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#    gunzip uniprot_sprot.fasta.gz
# 2. Make a blast database
#    makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot.fasta.db -dbtype prot
# 3. Run the blast
#    blastp -query [your fasta file] -db uniprot_sprot.fasta.db -outfmt 6 -max_target_seqs 1 -max_hsps 1 -evalue 1e-5 -out [your fasta file].blastp
# 4. Make the annotation file
#    python blast_annotation.py [your fasta file] [your blast output file]
# 5. Apply it to your phylogeny - open your phylogeny in FigTree, go to File>Import Annotations, and select your .annotation file.
#
# Usage: python blast_annotation.py [fasta_file] [blast_output]

import sys

# 1. load in fasta file
fasta = open(sys.argv[1],'r').read()
seq_d = {}
for line in fasta.split('>')[1:]:
    seq_d[line.split('\n')[0].strip().strip('>')] = []

# 2. load in blast file
blast = open(sys.argv[2],'r').readlines()
blast_d = {}
evalue_d = {}
for line in blast:
    if line.split('\t')[0].strip() in blast_d:
        pass
    else:
        blast_d[line.split('\t')[0].strip()] = line.split('\t')[1].split('_',1)[1].split('_',1)[1].split('OS=')[0].strip('_')
        evalue_d[line.split('\t')[0].strip()] = line.split('\t')[-2].strip()

# 3. Output annotation
out = open(sys.argv[1]+'.annotation','w')
out.write('seq\tprotein\tevalue\n')
for s in list(seq_d.keys()):
    try:
        out.write(s + '\t' + blast_d[s] + '\t' + evalue_d[s] + '\n')
    except:
        out.write(s + '\tno_annotation\tNA\n')
out.close()

