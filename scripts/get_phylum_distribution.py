#!/usr/bin/env python
import argparse
import pandas as pd
from tqdm import tqdm
from collections import Counter
from sys import argv

def parse_args(argv):
    desc = 'Collects metadata about antibiotic resistance genes.'
    parser = argparse.ArgumentParser(description=desc+'. ')
    parser.add_argument('--input', '-i', required=True,
                        help='Rarefied counts matrix produced by "rarefy_counts_matrix.py"')
    parser.add_argument('--cluster_dir', required=True,
                        help='Path to cluster directory.')
    parser.add_argument('--taxonomy_ncbi', required=True,
                        help='Table containing full lineage of bacterial genomes in NCBI.')
    parser.add_argument('--taxonomy_card', required=True,
                        help='Table containing full lineage of hosts carrying genes from the CARD database.')
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output file.')

    arguments = parser.parse_args()

    return arguments

def get_species(key, taxonomy_ncbi, taxonomy_card):
    if key.split('_')[0] == 'GCA':
        species = taxonomy_ncbi['_'.join(key.split('_')[:2])]
    else:
        species = taxonomy_card[key]

    return species

def lookup_host_species(arg, cluster_dir, taxonomy_ncbi, taxonomy_card):
    try:
        with open('/'.join([cluster_dir, arg.replace('(', '').replace(')', '').replace("'", "").replace('@',''), 'hidden.txt'])) as f:
            cluster_contents = [x.strip() for x in f.readlines()]

        species = [get_species(key, taxonomy_ncbi, taxonomy_card) for key in cluster_contents]

    except FileNotFoundError:
        species = 'NA'

    return species

def compare_lists(longer, shorter):
    shorter_set = set(shorter)

    return [1 if item in shorter_set else 0 for item in longer]

def read_dictionaries(rarefied_counts, cluster_dir, lineage_ncbi, lineage_card):
    taxonomy_ncbi = {}
    taxonomy_map = {}
    with open(lineage_ncbi, 'r') as f:
        for line in f:
            items = line.split('\t')
            taxonomy_ncbi[items[0]] = items[7].strip()
            taxonomy_map[items[7].strip()] = {'phylum': items[2], 'class': items[3], 'order': items[4], 'family': items[5]}


    taxonomy_card = {}
    with open(lineage_card, 'r') as f:
        for line in f:
            items = line.split('\t')
            taxonomy_card[items[0]] = ' '.join(items[1].split(' ')[:2])

    species_index = {}
    for arg in tqdm(rarefied_counts.columns):
        species_index[arg] = lookup_host_species(arg, cluster_dir, taxonomy_ncbi, taxonomy_card)

    return taxonomy_map, species_index

def get_host_phyla(lineage_ncbi):
    phyla = []
    with open(lineage_ncbi, 'r') as f:
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
    arguments = parse_args(argv)
    counts = pd.read_csv(arguments.input, sep='\t', index_col=0)
    taxonomy_map, species_index = read_dictionaries(counts, arguments.cluster_dir, arguments.taxonomy_ncbi, arguments.taxonomy_card)
    phyla = get_host_phyla(arguments.taxonomy_ncbi)

    arg_phylum_dist = count_host_phyla(counts, phyla, taxonomy_map, species_index)
    arg_phylum_dist.to_csv(arguments.output, index=True, header=True, sep='\t')

if __name__ == '__main__':
    main()