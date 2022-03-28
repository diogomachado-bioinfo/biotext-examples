#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
AUTHOR
Diogo de J. S. Machado

REQUIREMENTS
The testing of this script was done on Windows 10.
For this script to work it is necessary:
 - install the biotext library and its dependencies from PyPI. You can use the command "pip install biotext".
 - if you've never used the sweep library, you'll need to download the default projection matrix, just run: "python -c 'from sweep import down_proj_mat;down_proj_mat()'".
"""

import biotext as bt
import re
import pandas as pd

# Set this to true if you want to re-run the pubmed search,
# or keep it to False to use the previously downloaded set.
reRun_pubmedSearch = False
pubmedResult_fileName = 'pubmedResult_thioredoxin.tsv'
tree_fileName = 'tree_thioredoxin.nwk'

# Download or import pubmed result
if reRun_pubmedSearch:
    email = 'Your.Name.Here@example.org' # Email to inform pubmed who is doing the search.
    search_string = 'thioredoxin'
    
    pubmed_set = bt.pubSearch(search_string,email,batch_size=1000)
    pubmed_set.to_csv(pubmedResult_fileName,header=True,index=False,sep='\t')
else:
    pubmed_set = pd.read_csv(pubmedResult_fileName,sep='\t')

# Run aminocode
fasta = bt.aminocode.encodefile(pubmed_set['ti'],detailing='dp')

# Run SWeeP
mat = bt.fasta2vect(fasta)

# Create a string to identify the dendrogram leaves, using the PMIDs and paper titles.
# Non-alnumeric and non-space characters in titles are replaced with "_", to avoid conflict with the newick format syntax.
leaves = [str(r['pmid'])+' - '+re.sub('[^\w\s]','_',r['ti']) for i,r in pubmed_set.iterrows()]

# Run dendrogram (tree) creation in newick format
tree = bt.mat2tree(mat,leaves,method='complete')

# Save tree in a file
with open(tree_fileName,'w') as f:
    f.write(tree)
