#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
AUTHOR
Diogo de J. S. Machado

REQUIREMENTS
The testing of this script was done on Windows 10.
For this script to work it is necessary:
 - have the clustal omega available in the PATH of the operating system (download from <http://www.clustal.org/omega/#Download>);
 - install the used python libraries from PyPI (biotext, matplotlib and their dependencies);
 - if you've never used the sweep library, you'll need to download the default projection matrix, just run: "python -c 'from sweep import down_proj_mat;down_proj_mat()'".
"""

from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import string
from biotext.fastatools import (get_header,import_fasta,remove_pattern,
                                get_consensus,create_seqrecord_list,get_seq,
                                export_fasta,fasta_to_mat)
from biotext.aminocode import encode_list,decode_list
import os

# PARAMETERS
input_file_name = 'WP_011156533.1-02_27_2020.fasta'
aminocode_det = 'dp'
threshold = 0.015
plot_file_name = 'fastatext.png'
excel_file_name = 'fastatext.xlsx'
output_dir = 'fastaHeader2plot_result_dir'

try:
    os.mkdir(output_dir)
except:
   pass 

# read file
records = import_fasta(input_file_name)

# extract headers
descriptions = np.array(get_header(records))

# aminoencode
fastatxt = create_seqrecord_list(
    encode_list(descriptions,detail=aminocode_det))

# remove pattern
fastatxt = remove_pattern(fastatxt,['^(.*?)YS'])
     
# run sweep
mat_sweep = fasta_to_mat(fastatxt)

# PCA
pca = PCA(n_components=2)
pca.fit(mat_sweep.T)
pca_sweep=pca.components_.T * np.sqrt(pca.explained_variance_) # use loadings

# clustering
clus = AgglomerativeClustering(n_clusters = None, distance_threshold=threshold,
                               affinity='euclidean', linkage='ward')
clus.fit_predict(pca_sweep)
clus = clus.labels_

# define consensus
cons = []
for i in np.unique(clus):
    c, align = get_consensus(np.array(fastatxt, dtype=object)[clus==i])
    h = list(np.array(list(range(1,sum(clus==i)+1))).astype(str))
    export_fasta(
        create_seqrecord_list(align, header_list = h), output_dir+'/align_'+str(i)+
                  '.fasta') # save align file
    cons.append(c)
cons = create_seqrecord_list(cons)

# find and remove species name
cons = remove_pattern(cons,['YSYK\w+$']) # find and remove species name

# aminodecode
cons = get_seq(cons)
cons = decode_list(cons)

# create chart
fig, ax = plt.subplots()
for i in np.unique(clus):
    x = pca_sweep[clus==i,0]
    y = pca_sweep[clus==i,1]
    ax.scatter(x,y, label=cons[i]+' ('+str(len(x))+')', edgecolors='none')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.50, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('\nFactor loading 1', fontsize=12)
plt.ylabel('Factor loading 2\n', fontsize=12)
if isinstance(input_file_name, str):
    plt.title('Header text mining from file ' + input_file_name + '\n',
              fontsize=14)
else:
    plt.title('Header text mining\n', fontsize=14)
ax.grid(which='both')
if plot_file_name != None:
    plt.savefig(output_dir+'/'+plot_file_name, bbox_inches = 'tight', pad_inches = 0.1)

# save excel file to use in chart creating
if excel_file_name != None:
    alpha = list(string.ascii_uppercase)
    cons=np.array(cons)
    exFor = pd.DataFrame(columns=cons)
    for i in range(0,len(clus)):
        for ii in range(0,max(clus)+1):
            exFor.loc[i,cons[ii]] = ('=IF($A'+str(i+2)+'='+alpha[3+ii]+
                                     '$1,$C'+str(i+2)+',NA())')
    excel_df = pd.concat([pd.DataFrame(cons[clus],columns=['legend']),
                          pd.DataFrame(pca_sweep,columns=['a', 'b']),exFor],
                         axis=1)
    excel_df.to_excel(output_dir+'/'+excel_file_name, index=False)
