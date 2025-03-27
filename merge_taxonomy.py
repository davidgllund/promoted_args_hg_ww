import pandas as pd
import os
import fnmatch
from tqdm import tqdm
import re

def list_files(directory, filetype):
    """ Returns a list of files in 'directory' that match 'filetype' pattern. """
    return [file for file in os.listdir(directory) if fnmatch.fnmatch(file, filetype)]

def extract_taxid(column_name):
    match = re.search(r'\(taxid (\d+)\)', column_name)
    return match.group(1) if match else column_name

def merge_taxonomy_files(directory, output_file):
    """ Merges all *_taxonomy.csv files in the specified directory into one CSV file with consistent columns. """
    
    files = list_files(directory, "*.csv")
    if not files:
        print("No matching files found. Exiting.")
        return
    
    print(f"Found {len(files)} files. Merging...")

    dataframes = []
    
    for file in tqdm(files, desc="Processing files"):
        file_path = os.path.join(directory, file)
        try:
            df = pd.read_csv(file_path, index_col=None, sep=",", encoding="utf-8-sig", header=0)
            new_columns = [extract_taxid(x) for x in df.columns]
            df.columns = new_columns
            df = df.iloc[:,:len(df.columns)-2]
            df = df.T.groupby(df.columns).sum().T

            dataframes.append(df)
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
    
    # Merge all DataFrames
    merged_df = pd.concat(dataframes)

    merged_df.index = merged_df["TrueID"]
    merged_df.drop(columns=["TrueID"], inplace = True)

    print(f"Final merged DataFrame has {merged_df.shape[1]} columns and {merged_df.shape[0]} rows.")
    
    merged_df.to_csv(output_file, index=True, sep="\t")
    print(f"Merging complete! Output saved to {output_file}")

if __name__ == "__main__":
    merge_taxonomy_files("bacterial_taxonomy", "taxonomy_merged.tsv")