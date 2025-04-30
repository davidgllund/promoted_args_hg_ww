#!/usr/bin/env python
from sys import argv
import pandas as pd
import argparse
from tqdm import tqdm

def parse_args(argv):
    desc = 'Extract a subset of sequences from FASTA file.'
    parser = argparse.ArgumentParser(description=desc+'. ')
    parser.add_argument('--input', '-i', required=True,
                        help='Input file (BLAST output format 6).')
    parser.add_argument('--seq_lengths', required=True,
                        help='File the lengths of all query sequences.')
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output file.')
    
    arguments = parser.parse_args()
    
    return arguments

def pre_filter_blast_results(blastout):
    blastout.loc[:,12] = blastout.iloc[:,7]-blastout.iloc[:,6]
    blastout = blastout.loc[blastout.iloc[:,2] >= 90, :]
    blastout.loc[:,0] = blastout.loc[:,0].astype(str)

    return blastout

def find_best_blast_hit(blastout, subject_lengths):
    tqdm.pandas()
    
    blastout = blastout.merge(subject_lengths, left_on=1, right_index=True, how='left')
    blastout.columns = list(blastout.columns[:-1]) + ['subject_length']

    blastout['cover_threshold'] = blastout['subject_length'] * 0.7

    blastout = blastout[blastout.iloc[:, 12] >= blastout['cover_threshold']]

    best_hits = (
        blastout.groupby(blastout.iloc[:, 0], group_keys=False)
        .progress_apply(lambda g: g.loc[g[2].idxmax()])
    )

    return best_hits

def main():
    arguments = parse_args(argv)
    
    print("Reading BLAST output file...")
    blast_results = pd.read_csv(arguments.input, delimiter='\t', header=None)

    print("Pre-filtering BLAST output file...")
    blast_results = pre_filter_blast_results(blast_results)

    print("Reading sequence lengths file...")
    seq_lengths = pd.read_csv(arguments.seq_lengths, sep="\t", header=None, index_col=0)

    print("Searching for best BLAST hits...")
    best_hits = find_best_blast_hit(blast_results, seq_lengths)

    best_hits.to_csv(arguments.output, sep="\t", index=False)

if __name__ == '__main__':
    main()
