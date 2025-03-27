#!/usr/bin/env python
import pandas as pd
from tqdm import tqdm
import numpy as np
from collections import Counter

def get_species(key, taxonomy_ncbi, taxonomy_card):
    if key.split('_')[0] == 'GCA':
        species = taxonomy_ncbi['_'.join(key.split('_')[:2])]
    else:
        species = taxonomy_card[key]

    return species

def lookup_host_species(arg, taxonomy_ncbi, taxonomy_card):
    try:
        with open('/'.join(['/storage/dlund/Michaela_reloaded/data/clusters', arg.replace('(', '').replace(')', '').replace("'", "").replace('@',''), 'hidden.txt'])) as f:
            cluster_contents = [x.strip() for x in f.readlines()]

        species = [get_species(key, taxonomy_ncbi, taxonomy_card) for key in cluster_contents]

    except FileNotFoundError:
        species = 'NA'

    return species

def read_dictionaries(counts):
    taxonomy_ncbi = {}
    with open('/home/dlund/index_files/genome_full_lineage.tsv', 'r') as f:
        for line in f:
            items = line.split('\t')
            taxonomy_ncbi[items[0]] = items[7].strip()

    taxonomy_card = {}
    with open('/home/dlund/index_files/card_taxonomy.txt', 'r') as f:
        for line in f:
            items = line.split('\t')
            taxonomy_card[items[0]] = ' '.join(items[1].split(' ')[:2])

    species_index = {}
    for arg in tqdm(counts.columns):
        species_index[arg] = lookup_host_species(arg, taxonomy_ncbi, taxonomy_card)

    return species_index

def count_pathogens(species_list, pathogen_list):
    hits = []

    species_count = Counter(species_list)

    if len(species_list) >= 3:
            threshold = 3
    else:
        threshold = len(species_list)

    host_species = [item for item, count in species_count.items() if count >= threshold and item != "NA"]
    
    for sp in pathogen_list:
        if sp in host_species:
            hits.append(1)
        else:
            hits.append(0)

    return hits

def main():
    counts = pd.read_csv('rarefied_counts_updated.tsv', sep='\t', index_col=0)
    species_index = read_dictionaries(counts)

    with open('/home/dlund/index_files/pathogen_list.txt', 'r') as f:
        pathogen_index = [x.strip() for x in f.readlines()]

    pathogens = [count_pathogens(species_index[arg], pathogen_index) for arg in tqdm(counts.columns)]
    
    df = pd.DataFrame(pathogens, columns=pathogen_index, index=counts.columns).transpose()
    df.to_csv("arg_pathogen_distr.tsv", sep="\t", header=True)

if __name__ == '__main__':
    main()