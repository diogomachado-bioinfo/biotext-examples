#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
AUTHOR
Diogo de J. S. Machado

REQUIREMENTS
The testing of this script was done on Windows 10.
For this script to work it is necessary:
 - have the clustal omega available in the PATH of the operating system (downlaod from <http://www.clustal.org/omega/#Download>);
 - install the used python libraries from PyPI (biotext, sweep, sklearn, matplotlib and pandas).
 
RECOMMENDATION
Files will be generated with the execution of this script, so it is recommended that it be executed inside an empty folder, for organizational reasons.
"""

from sweep import fas2sweep
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import string
from biotext.fastatools import getHeader, fastaread, removePattern, getCons, list2fasta, getSeq
from biotext.aminocode import encodefile, decodetext

# PARAMETERS
input_file_name = 'WP_011156533.1-02_27_2020.fasta'
aminocode_det = 'dp'
threshold = 0.015
plot_file_name = 'fastatext.png'
excel_file_name = 'fastatext.xlsx'

# read file
records = fastaread(input_file_name)

# extract headers
descriptions = np.array(getHeader(records))

# aminoencode
fastatxt = encodefile(descriptions,detailing=aminocode_det)

# remove pattern
fastatxt = removePattern(fastatxt,['^(.*?)YS'])
     
# sweep
mat_sweep = fas2sweep(fastatxt)

# PCA
pca = PCA(n_components=2)
pca.fit(mat_sweep.T)
pca_sweep=pca.components_.T * np.sqrt(pca.explained_variance_) # use loadings

# clustering
clus = AgglomerativeClustering(n_clusters = None, distance_threshold=threshold, affinity='euclidean', linkage='ward')
clus.fit_predict(pca_sweep)
clus = clus.labels_

# define consensus
cons = []
for i in np.unique(clus):
    c, align = getCons(np.array(fastatxt)[clus==i])
    h = list(np.array(list(range(1,sum(clus==i)+1))).astype(str))
    fastawrite (list2fasta(align, header = h), 'align_'+str(i)+'.fasta') # save align file
    cons.append(c)
cons = list2fasta(cons)

# find and remove species name
cons = removePattern(cons,['YSYK\w+$']) # find and remove species name

# aminodecode
cons = getSeq(cons)
dt_vect = np.vectorize(lambda x: decodetext(x,detailing=aminocode_det))
cons = dt_vect(cons)

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
    plt.title('Header text mining from file ' + input_file_name + '\n', fontsize=14)
else:
    plt.title('Header text mining\n', fontsize=14)
ax.grid(which='both')
if plot_file_name != None:
    plt.savefig(plot_file_name, bbox_inches = 'tight', pad_inches = 0.1)

# save excel file to use in chart creating
if excel_file_name != None:
    alpha = list(string.ascii_uppercase)
    cons=np.array(cons)
    exFor = pd.DataFrame(columns=cons)
    for i in range(0,len(clus)):
        for ii in range(0,max(clus)+1):
            exFor.loc[i,cons[ii]] = '=IF($A'+str(i+2)+'='+alpha[3+ii]+'$1,$C'+str(i+2)+',NA())'
    excel_df = pd.concat([pd.DataFrame(cons[clus],columns=['legend']),pd.DataFrame(pca_sweep,columns=['a', 'b']),exFor],axis=1)
    excel_df.to_excel(excel_file_name, index=False)
