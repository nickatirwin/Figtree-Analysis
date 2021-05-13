# ANNOTATE A PHYLOGENY WITH TAXONOMY INFORMATION
#
# Create an annotation file for FigTree from taxa ids in the headers of a fasta
# 
# NOTE: fasta file headers must be in the form "NCBITaxaID.ProteinID" (e.g. 9606.Q53XC5)
#
'''
# install ete3
conda install -c bioconda ete3

# download and update the NCBI Taxonomy database (run in python)
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

# run the script
python ncbi_taxid_to_annotation.py ProteinA.fasta
'''

from ete3 import NCBITaxa
ncbi = NCBITaxa()
import sys
import subprocess

#load in the fasta file
lines = open(sys.argv[1], 'r').readlines()
output = open(sys.argv[1] + '.taxonomy_annotation', 'w')
otu = []

#extract the headers (OTUs) and taxa ids
taxaids = []
for i in lines:
    if i.startswith('>'):
        taxaids.append(i.split('.')[0].strip('>').strip())
        otu.append(i.strip('>').strip())

#get the names corresponding to each of the taxa ids using grep
taxon = []
domain = []
for i in taxaids:
    try:
        taxon.append(ncbi.get_taxid_translator([int(i)])[int(i)])
        if 2759 in ncbi.get_lineage(int(i)):
            domain.append('Eukaryota')
        elif 10239 in ncbi.get_lineage(int(i)):
            domain.append('Virus')
        elif 2 in ncbi.get_lineage(int(i)):
            domain.append('Bacteria')
        else:
            domain.append('Archaea')
        supergroup.append(str(ncbi.get_taxid_translator([ncbi.get_lineage(int(i))[3]])[ncbi.get_lineage(int(i))[3]]))
    except:
        taxon.append(i)
        domain.append(i)
        supergroup.append(i)
#output the results
output.write('otu\ttaxaid\tspecies\tsupergroup\tdomain\n')
n = 0
while n < len(taxaids):
    output.write(otu[n].strip('\n') + '\t' + taxaids[n] + '\t' + taxon[n].strip('\n') + '\t' + supergroup[n] + '\t' +  domain[n].strip() + '\n')
    n += 1

output.close()
