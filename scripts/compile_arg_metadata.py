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
    parser.add_argument('--arg_index', required=True,
                        help='Table containing gene class and antibiotic for the included ARGs.')
    parser.add_argument('--blastout', required=True,
                        help='BLAST output file showing similiarity to genes in the ResFinder database.')
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

def latent_or_established(id, index, flag):
    if id.split('_')[0] != 'GCA':
        cat = 'Established'
        label = id

    else:
        if id in index.keys():
            cat = 'Established'
            label = index[id]
        else:
            cat = 'Latent'
            label = id

    if flag == 'category':
        return cat
    
    elif flag == 'label':
        return label
    
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

def get_taxonomy_score(species_list, taxonomy_map):
    try:
        phylum_list = Counter([taxonomy_map[sp]['phylum'] for sp in species_list if sp in taxonomy_map.keys()])
        class_list = Counter([taxonomy_map[sp]['class'] for sp in species_list if sp in taxonomy_map.keys()])
        order_list = Counter([taxonomy_map[sp]['order'] for sp in species_list if sp in taxonomy_map.keys()])
        family_list = Counter([taxonomy_map[sp]['family'] for sp in species_list if sp in taxonomy_map.keys()])
        genus_list = Counter([x.split(' ')[0] for x in species_list])

        if len(species_list) >= 3:
            threshold = 3
        else:
            threshold = len(species_list)

        host_phyla = [item for item, count in phylum_list.items() if count >= threshold and item != "NA"]
        host_classes = [item for item, count in class_list.items() if count >= threshold and item != "NA"]
        host_orders = [item for item, count in order_list.items() if count >= threshold and item != "NA"]
        host_families = [item for item, count in family_list.items() if count >= threshold and item != "NA"]
        host_genera = [item for item, count in genus_list.items() if count >= threshold and item != "NA"]

        if len(host_phyla) > 1:
            score = 6
        elif len(host_phyla) == 1 and len(host_classes) > 1:
            score = 5
        elif len(host_classes) == 1 and len(host_orders) > 1:
            score = 4
        elif len(host_orders) == 1 and len(host_families) > 1:
            score = 3
        elif len(host_families) == 1 and len(host_genera) > 1:
            score = 2
        elif len(host_genera) == 1 and len(species_list) > 1:
            score = 1
        elif len(species_list) == 1:
            score = 0
        else:
            score = 'NA'

    except KeyError:
        score = 'NA'
            
    return score

def count_pathogens(species_list, pathogen_list):
    n = 0
    species_count = Counter(species_list)

    if len(species_list) >= 3:
            threshold = 3
    else:
        threshold = len(species_list)

    host_species = [item for item, count in species_count.items() if count >= threshold and item != "NA"]
    
    for sp in host_species:
        if sp in pathogen_list:
            n += 1

    return n

def count_phyla(species_list, taxonomy_map):
    try:
        phylum_list = Counter([taxonomy_map[sp]['phylum'] for sp in species_list if sp in taxonomy_map.keys()])

        if len(species_list) >= 3:
            threshold = 3
        else:
            threshold = len(species_list)

        host_phyla = [item for item, count in phylum_list.items() if count >= threshold and item != "NA"]
        n = len(set(host_phyla))
        
    except KeyError:
        n = 1
    
    return n

def read_dictionaries(rarefied_counts, index, blastout, cluster_dir, lineage_ncbi, lineage_card):
    arg_index = {}
    with open(index, 'r') as f:
        for line in f:
            items = line.split('\t')
            arg_index[items[0]] = items[1]

    established_index = {}
    with open(blastout, 'r') as f:
        for line in f:
            items = line.split('\t')
            established_index[items[0]] = '_'.join(items[1].split('_')[:-3])

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
    for arg in tqdm(rarefied_counts.columns, desc='Collecting taxonomic information.'):
        species_index[arg] = lookup_host_species(arg, cluster_dir, taxonomy_ncbi, taxonomy_card)

    return arg_index, established_index, taxonomy_map, species_index

def main():
    arguments = parse_args(argv)
    counts = pd.read_csv(arguments.input, sep='\t', index_col=0)
    
    arg_index, established_index, taxonomy_map, species_index = read_dictionaries(counts, arguments.cluster_dir, arguments.arg_index, arguments.blastout, arguments.taxonomy_ncbi, arguments.taxonomy_card)
    
    with open(arguments.pathogens, 'r') as f:
        pathogen_index = [x.strip() for x in f.readlines()]

    counts.loc['arg_class'] = [arg_index[key] for key in counts.columns]
    counts.loc['category'] = [latent_or_established(key, established_index, 'category') for key in counts.columns]
    counts.loc['label'] = [latent_or_established(key, established_index, 'label') for key in counts.columns]
    counts.loc['taxonomy_score'] = [get_taxonomy_score(species_index[arg], taxonomy_map) for arg in tqdm(counts.columns, desc='Estaimating taxonomic range.')]
    counts.loc['n_pathogens'] = [count_pathogens(species_index[arg], pathogen_index) for arg in tqdm(counts.columns, desc='Counting pathogenic bacterial hosts.')]
    counts.loc['n_phyla'] = [count_phyla(species_index[arg], taxonomy_map) for arg in tqdm(counts.columns, desc='Counting host phyla.')]
    counts.loc['cluster_size'] = [len(species_index[arg]) for arg in tqdm(counts.columns, desc='Calculating cluster size.')]

    metadata_rows = ['arg_class', 'category', 'label', 'cluster_size', 'taxonomy_score', 'n_phyla', 'n_pathogens']
    metadata = counts.loc[metadata_rows]
    counts = counts.loc[list(set(counts.index) - set(metadata_rows))]

    metadata.to_csv(arguments.output, sep = '\t', header=True, index=True)

if __name__ == '__main__':
    main()