#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from sys import argv

def parse_args(argv):
    desc = 'Rarefies counts matrix.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--input', '-i', required=True,
                        help='Counts matrix produced by "compile_counts_matrix.py".')
    parser.add_argument('--n_reads', required=True,
                        help='List of the total number of reads included in each sample.')
    parser.add_argument('--subset', required=True,  help='List of ARGs to analyze.'),
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output file.')

    arguments = parser.parse_args()

    return arguments

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

def subset_counts_matrix(counts, subset):
    counts_compiled = pd.read_csv(counts, sep='\t', index_col=0)

    with open(subset, 'r') as f:
        subset = [x.strip() for x in f.readlines()]

    counts_compiled = counts_compiled[subset]

    return counts_compiled

def add_non_arg_reads(counts, total_reads):
    read_counts = {}
    with open(total_reads, 'r') as f:
        for line in f:
            items = line.split('\t')
            read_counts[items[0]] = int(items[1].strip())

    non_arg = [get_non_matched_reads(key, counts, read_counts) for key in tqdm(counts.index)]
    counts['non-ARG'] = non_arg

    return counts

def main():
    arguments = parse_args(argv)

    counts_compiled = subset_counts_matrix(arguments.input, arguments.subset)
    counts_compiled = add_non_arg_reads(counts_compiled, arguments.n_reads)

    print('Rarefying data')
    rarefied = [rarefy(key, counts_compiled, 5000000) for key in tqdm(counts_compiled.index)]
    rarefied_counts = pd.DataFrame(data=rarefied, columns=counts_compiled.columns, index=counts_compiled.index)
    rarefied_counts = rarefied_counts.drop('non-ARG', axis=1)
    rarefied_counts.to_csv(arguments.output, sep='\t', header=True, index=True)

if __name__ == '__main__':
    main()