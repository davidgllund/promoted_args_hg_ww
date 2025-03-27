#!/usr/bin/env python
import pandas as pd
from tqdm import tqdm
from itertools import product
from collections import Counter
from sys import argv
import math
import argparse

def parse_args(argv):
    desc = 'Compile info about co-localized mobile genetic elements'
    copyright = 'Copyright (c) David Lund 2024.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--input', '-i', required=True,
                        help='Fasta file containing genes to analyze.')
    parser.add_argument('--genome_kmers', required=True,
                        help='File containing the kmer distributions of relevant genomes.')
    parser.add_argument('--pathogens', required=True,
                        help='File containing pathogenic species and representative NCBI Assembly IDs.')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file name.')
    
    arguments = parser.parse_args()
    
    return arguments

def read_genome_kmers(file_path):
    kmer_distributions = {}

    with open(file_path, 'r') as f:
        total_lines = sum(1 for _ in f)

    with open(file_path, 'r') as f:
        next(f)

        for line in tqdm(f, total=total_lines - 1, desc="Loading k-mer distributions"):
            items = line.split('\t')
            kmer_distributions[items[0]] = [float(x) for x in items[1:]]

    return kmer_distributions

def fasta_reader(filename):
    fasta_df = pd.read_csv(filename, sep='>', lineterminator='>', header=None)
    fasta_df[['Accession', 'Sequence']] = fasta_df[0].str.split('\n', n=1, expand=True)
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df['Sequence'] = fasta_df['Sequence'].replace('\n', '', regex=True)
    fasta_df.index = list(fasta_df['Accession'])
    
    return fasta_df

def get_kmers(seq, k):
    kmers = []
    
    for i in range(len(seq)-(k-1)):
        kmers.append(str(seq[i:i+k]))
        
    return kmers

def generate_possible_kmers(k):
    comb = list(product('ACGT', repeat=k))
    possible_kmers = []

    for i in range(len(comb)):
        possible_kmers.append(''.join(comb[i]))

    return possible_kmers

def get_kmer_distribution(kmers, possible_kmers):
    kmers_counted = Counter(kmers)
    dist_list = []

    for kmer in possible_kmers:
        dist_list.append(kmers_counted[kmer]/len(kmers))

    distribution = pd.DataFrame({'fraction': dist_list}, index=possible_kmers)

    return list(distribution['fraction'])

def main():
    arguments = parse_args(argv)

    args = fasta_reader(arguments.input)
    kmer_distributions = read_genome_kmers(arguments.genome_kmers)
    pathogens = pd.read_csv(arguments.pathogens, sep="\t", index_col=0, header=None)
    
    all_kmers = generate_possible_kmers(5)

    gene_kmers = {}
    for arg in tqdm(args.index, desc="Calculating gene k-mer distributions"):
        kmers = get_kmers(args.loc[arg, 'Sequence'], 5)
        gene_kmers[arg] = get_kmer_distribution(kmers, all_kmers)

    results = pd.DataFrame(index=pathogens.index, columns=gene_kmers.keys())
    for sp in tqdm(results.index, desc="Computing k-mer distances"):
        results.loc[sp] = [math.dist(kmer_distributions[list(pathogens.loc[sp])[0]], gene_kmers[x]) for x in results.columns]

    results.to_csv(arguments.output, sep='\t', header=True, index=True)

if __name__ == '__main__':
    main()
