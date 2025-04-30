#!/usr/bin/env python
import pandas as pd
import numpy as np
from tqdm import tqdm

def get_non_matched_reads(sample, arg_counts, total_counts):
    non_matched = total_counts[sample] - sum(arg_counts.loc[sample,:])
    
    return non_matched

def rarefy(id, counts, depth):
    counts_row = counts.loc[id]
    reads = np.repeat(counts_row.index, counts_row.values)
    
    reads_sample = np.random.choice(reads, size=depth, replace=False)
    
    counts_sample = pd.Series(reads_sample).value_counts()
    
    rarefied_counts = pd.Series(0, index=counts.columns)
    rarefied_counts.update(counts_sample)
    
    return rarefied_counts.values.tolist()

def subset_counts_matrix():
    counts_compiled = pd.read_csv('counts_compiled.tsv', sep='\t', index_col=0)
    with open('../relevant_subset.txt', 'r') as f:
        subset = [x.strip() for x in f.readlines()]

    counts_compiled = counts_compiled[subset]

    return counts_compiled

def add_non_arg_reads(counts):
    read_counts = {}
    with open('../QC_read_count_ebi.txt', 'r') as f:
        for line in f:
            items = line.split('\t')
            if items[1].split('_')[1] == '1':
                read_counts[items[1].split('_')[0]] = int(items[2].strip())
            else:
                continue

    non_arg = [get_non_matched_reads(key, counts, read_counts) for key in tqdm(counts.index)]
    counts['non-ARG'] = non_arg

    return counts

def main():
    print("Compiling counts matrix")
    counts_compiled = subset_counts_matrix()
    counts_compiled = add_non_arg_reads(counts_compiled)

    print('Rarefying data')
    rarefied = [rarefy(key, counts_compiled, 5000000) for key in tqdm(counts_compiled.index)]
    rarefied_counts = pd.DataFrame(data=rarefied, columns=counts_compiled.columns, index=counts_compiled.index)
    rarefied_counts = rarefied_counts.drop('non-ARG', axis=1)
    rarefied_counts.to_csv('rarefied_counts.tsv', sep='\t', header=True, index=True)

if __name__ == '__main__':
    main()