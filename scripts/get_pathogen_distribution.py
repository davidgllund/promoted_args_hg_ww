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
    parser.add_argument('--pathogens', required=True,
                        help='List of pathogenic bacterial species.')
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

def read_dictionaries(rarefied_counts, cluster_dir, lineage_ncbi, lineage_card):
    taxonomy_ncbi = {}
    with open(lineage_ncbi, 'r') as f:
        for line in f:
            items = line.split('\t')
            taxonomy_ncbi[items[0]] = items[7].strip()

    taxonomy_card = {}
    with open(lineage_card, 'r') as f:
        for line in f:
            items = line.split('\t')
            taxonomy_card[items[0]] = ' '.join(items[1].split(' ')[:2])

    species_index = {}
    for arg in tqdm(rarefied_counts.columns):
        species_index[arg] = lookup_host_species(arg, cluster_dir, taxonomy_ncbi, taxonomy_card)

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
    arguments = parse_args(argv)
    counts = pd.read_csv(arguments.input, sep='\t', index_col=0)
    species_index = read_dictionaries(counts, arguments.cluster_dir, arguments.taxonomy_ncbi, arguments.taxonomy_card)

    with open(arguments.pathogens, 'r') as f:
        pathogen_index = [x.strip() for x in f.readlines()]

    pathogens = [count_pathogens(species_index[arg], pathogen_index) for arg in tqdm(counts.columns)]
    
    df = pd.DataFrame(pathogens, columns=pathogen_index, index=counts.columns).transpose()
    df.to_csv(arguments.output, sep="\t", header=True)

if __name__ == '__main__':
    main()