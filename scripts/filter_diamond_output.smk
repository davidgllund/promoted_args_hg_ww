#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
# IMPORT PACKAGES AND DEFINE FUNCTIONS
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import random

directories, files = glob_wildcards('example_data/{dir}/{file}.tsv')

def prefilter_data(filename):
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen']

    try:
        df = pd.read_csv(filename, sep='\t', header=None)
        df.columns = cols
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns = cols)

    df = df.loc[df.loc[:, 'length'] >= 20, :]
    df = df.loc[df.loc[:, 'length'] / df.loc[:, 'qlen'] >= 0.20, :]

    return df

def select_best_hit(query_id, dataframe):
    subset = dataframe.loc[dataframe.iloc[:,0] == query_id.iloc[0], :]
    best_id = subset.loc[subset.loc[:, 'pident'] == max(subset.loc[:,'pident']), :]
    best_length = best_id.loc[best_id.loc[:, 'length'] == max(best_id.loc[:,'length']), :]
    best_eval = best_length.loc[best_length.loc[:, 'evalue'] == max(best_length.loc[:,'evalue']), :]

    if len(best_eval.iloc[:,0]) > 1:
        best_random = best_eval.iloc[random.sample(range(len(best_eval)), 1),:]

        return best_random.iloc[0,:]
    
    else:
        return best_eval.iloc[0,:]

rule all:
    input:
        expand('example_data/{dir}/filtered/{file}_filtered.txt',
        zip, dir=directories, file=files)

rule filter:
    input:
        'example_data/{dir}/{file}.tsv'
    output:
        'example_data/{dir}/filtered/{file}_filtered.txt'
    run:
        df = prefilter_data(input[0])
        query_ids = pd.DataFrame({'id': df.iloc[:,0].unique()})

        best_hits = query_ids.apply(select_best_hit, dataframe=df, axis=1)
        best_hits.to_csv(output[0], sep='\t', index=False)

    