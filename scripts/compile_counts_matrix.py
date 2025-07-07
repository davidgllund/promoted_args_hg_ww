#!/usr/bin/env python
import pandas as pd
import argparse
from sys import argv
from tqdm import tqdm

def parse_args(argv):
    desc = 'Creates counts matrix from DIAMOND output files.'
    parser = argparse.ArgumentParser(description=desc+'. ')
    parser.add_argument('--files', required=True,
                        help='List of paths to DIAMOND output files.')
    parser.add_argument('--ids', required=True,
                        help='List of genes to include.')
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output file.')

    arguments = parser.parse_args()

    return arguments

def compile_counts(paths, ids):
    counts_total = []
    new_index = []

    for i in tqdm(range(len(paths)), desc = 'Compiling counts matrix'):
        df = pd.read_csv(paths.iloc[i,0], sep='\t')

        if (len(df) == 0):
            continue
        else:
            df_counted = pd.DataFrame(df['sseqid'].value_counts())

            sample_counts = []
            for id in ids:
                if id in list(df_counted.index):
                    sample_counts.append(df_counted.loc[id, 'count'])
                else:
                    sample_counts.append(0)

            counts_total.append(sample_counts)
            new_index.append(list(paths.index)[i])

    counts_compiled = pd.DataFrame(data=counts_total)
    counts_compiled.columns = ids
    counts_compiled.index = new_index

    return counts_compiled

def main():
    arguments = parse_args(argv)
    paths = pd.read_csv(arguments.files, sep='\t', index_col=0, header=None)

    with open(arguments.ids, 'r') as f:
        ids = [x.strip() for x  in f.readlines()]

    counts_compiled = compile_counts(paths, ids)    
    counts_compiled.to_csv(arguments.output, sep='\t', header=True, index=True)

if __name__ == '__main__':
    main()