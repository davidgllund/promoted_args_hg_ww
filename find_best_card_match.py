#!/usr/bin/env python
#-------------------------------------------------------------------------------
# DESCRIPTION
#-------------------------------------------------------------------------------
# This script analyzes BLAST results to find the best match between OTUs and 
# bacterial genomes.
#
# INPUT ARGUMENTS:
# 1. BLAST output file in standard format 6.
# 2. File containing lengths of subject sequences. 
# 3. Name of output file.
#
#-------------------------------------------------------------------------------
# IMPORT PACKAGES
#-------------------------------------------------------------------------------
import sys
import pandas as pd
import numpy as np
from itertools import islice
import progressbar
import time

#-------------------------------------------------------------------------------
# PREPARE FILES
#-------------------------------------------------------------------------------
print("Reading BLAST output file")
blast_results = pd.read_csv("blastout_card.txt", delimiter='\t', header=None)

print("Pre-filtering BLAST output file")
blast_results.loc[:,12] = blast_results.iloc[:,7]-blast_results.iloc[:,6]
blast_results = blast_results.loc[blast_results.iloc[:,2] >= 90, :]
blast_results.loc[:,0] = blast_results.loc[:,0].astype(str)

seq_lengths = pd.read_csv(sys.argv[2], sep="\t", header=None, index_col=0)
#-------------------------------------------------------------------------------
# ANALYSIS
#-------------------------------------------------------------------------------
# For each subject, find the hit with highest sequence identity
query_ids = blast_results.iloc[:,0].unique()
df = pd.DataFrame()

print("Searching for best BLAST hits")

widgets = [' [',
            progressbar.Timer(format= 'elapsed time: %(elapsed)s'),
            '] ',
            progressbar.Bar('*'),' (',
            progressbar.ETA(), ') ',
            ]

bar = progressbar.ProgressBar(max_value=len(subject_ids), widgets=widgets).start()

for i in range(len(query_ids)):
    subset = blast_results.loc[blast_results.iloc[:,0] == query_ids[i], :]
    length = seq_lengths.loc[subset.iloc[0,1]]
    cover_threshold = length.iloc[0] * 0.7
    subsubset = subset.loc[pd.Series(subset.iloc[:,12] >= float(cover_threshold)), :] 
    
    if len(subsubset.iloc[:,0]) == 0:

        time.sleep(0.1)
        bar.update(i)
        
        continue
    
    else:
        df = df.append(subsubset.iloc[np.argmax(subsubset.iloc[:,2]), :])

        time.sleep(0.1)
        bar.update(i)
        
# Write results to file
df.to_csv(sys.argv[3])

#-------------------------------------------------------------------------------
