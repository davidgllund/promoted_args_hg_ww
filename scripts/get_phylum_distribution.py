#!/usr/bin/env python
import pandas as pd
from tqdm import tqdm
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

def compare_lists(longer, shorter):
    shorter_set = set(shorter)

    return [1 if item in shorter_set else 0 for item in longer]

def read_dictionaries(counts):
    taxonomy_ncbi = {}
    taxonomy_map = {}
    with open('/home/dlund/index_files/genome_full_lineage.tsv', 'r') as f:
        for line in f:
            items = line.split('\t')
            taxonomy_ncbi[items[0]] = items[7].strip()
            taxonomy_map[items[7].strip()] = {'phylum': items[2], 'class': items[3], 'order': items[4], 'family': items[5]}

    taxonomy_card = {}
    with open('/home/dlund/index_files/card_taxonomy.txt', 'r') as f:
        for line in f:
            items = line.split('\t')
            taxonomy_card[items[0]] = ' '.join(items[1].split(' ')[:2])

    species_index = {}
    for arg in tqdm(counts.columns):
        species_index[arg] = lookup_host_species(arg, taxonomy_ncbi, taxonomy_card)

    return taxonomy_map, species_index

def get_host_phyla():
    phyla = []
    with open('/home/dlund/index_files/genome_full_lineage.tsv', 'r') as f:
        for line in f:
            items = line.split('\t')
            phyla.append(items[2])

    phyla = list(set(phyla))
    phyla.sort()

    return phyla

def count_host_phyla(counts, phyla, taxonomy_map, species_index):
    arg_phylum_dist = pd.DataFrame(data = 0, columns = counts.columns, index = phyla)

    for arg in tqdm(arg_phylum_dist.columns):
        phylum_list = [taxonomy_map[x]['phylum'] for x in species_index[arg] if x in taxonomy_map.keys()]

        if len(phylum_list) >= 3:
            threshold = 3
        else:
            threshold = len(phylum_list)

        phylum_counts = Counter(phylum_list)
        filtered_phyla = [item for item, count in phylum_counts.items() if count >= threshold and item != "NA"]
        arg_phylum_dist[arg] = compare_lists(arg_phylum_dist.index, filtered_phyla)

    return arg_phylum_dist

def main():
    counts = pd.read_csv('rarefied_counts_updated.tsv', sep='\t', index_col=0)
    taxonomy_map, species_index = read_dictionaries(counts)
    phyla = get_host_phyla()

    arg_phylum_dist = count_host_phyla(counts, phyla, taxonomy_map, species_index)
    arg_phylum_dist.to_csv('arg_phylum_distr.tsv', index=True, header=True, sep='\t')

if __name__ == '__main__':
    main()