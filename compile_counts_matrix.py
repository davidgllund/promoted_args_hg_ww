#!/usr/bin/env python
import pandas as pd

def compile_counts(paths, ids):
    counts_total = []
    new_index = []

    for i in range(len(paths)):
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
            print(i)

    counts_compiled = pd.DataFrame(data=counts_total)
    counts_compiled.columns = ids
    counts_compiled.index = new_index

    return counts_compiled

def main():
    paths = pd.read_csv('sample_paths.tsv', sep='\t', index_col=0, header=None)
    
    with open('../arg_ids.txt', 'r') as f:
        ids = [x.strip() for x  in f.readlines()]

    counts_compiled = compile_counts(paths, ids)    
    counts_compiled.to_csv('counts_compiled.tsv', sep='\t', header=True, index=True)

if __name__ == '__main__':
    main()