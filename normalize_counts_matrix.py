#!/usr/bin/env python
import pandas as pd
import numpy as np
from tqdm import tqdm

def normalize(counts_matrix, counts_index, sample):
    div = (counts_matrix.loc[sample]+1) / counts_index[sample]
    normalized = np.log(10**6 * div)
    
    return normalized

def rarefy_bacterial_counts():
    bacterial_counts = {}
    with open('bacterial_counts.tsv', 'r') as f:
        for line in f:
            items = line.split('\t')
            bacterial_counts[items[0].split('.')[0]] = int(items[1].strip())

    non_bacterial_count = {}
    with open('qc_read_count.tsv', 'r') as f:
        for line in f:
            items = line.split('\t')
            if items[0] in bacterial_counts.keys():
                non_bacterial_count[items[0]] = int(items[1]) - bacterial_counts[items[0]]

    rarefied_bacteria = {}
    for key in tqdm(bacterial_counts.keys()):
        if key in non_bacterial_count.keys():
            bacteria_count = bacterial_counts[key]
            other_count = non_bacterial_count[key]
            total_count = bacteria_count + other_count

            if total_count >= 5000000:
                chosen_indices = np.random.choice(total_count, size=5000000, replace=False)
                rarefied_bacteria[key] = np.sum(chosen_indices < bacteria_count)

    return rarefied_bacteria

def main():
    counts = pd.read_csv('rarefied_counts.tsv', sep='\t', index_col=0)

    with open('../relevant_subset.txt', 'r') as f:
        subset = [x.strip() for x in f.readlines()]

    counts = counts[subset]

    rarefied_bacteria = rarefy_bacterial_counts()
    normalized = [normalize(counts, rarefied_bacteria, sample) for sample in tqdm(counts.index) if sample in rarefied_bacteria.keys()]

    included_samples = [x for x in counts.index if x in rarefied_bacteria.keys()]

    normalized_counts = pd.DataFrame(data=normalized, columns=counts.columns, index=included_samples)
    normalized_counts[counts == 0] = 0
    normalized_counts.to_csv("normalized_counts_human.tsv", sep="\t", header=True, index=True)

if __name__ == '__main__':
    main()
