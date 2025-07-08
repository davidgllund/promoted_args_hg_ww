#!/usr/bin/bash

# Unpack and restructure data
bash scripts/unpack_data.sh

# Filter DIAMOND output files
snakemake -s scripts/filter_diamond_output.smk --cores 2 all
mkdir tmp

# Compile and rarefy counts matrix (human gut samples)
find example_data/mapped_reads_human_gut/filtered -type f > tmp/paths.txt
cat tmp/paths.txt | rev |cut -d '/' -f 1 | rev | cut -d '_' -f 1 > tmp/ids.txt
paste tmp/ids.txt tmp/paths.txt > tmp/sample_paths_hg.txt
rm tmp/paths.txt tmp/ids.txt

python scripts/compile_counts_matrix.py --files tmp/sample_paths_hg.txt --ids auxiliary_files/arg_ids.txt --output tmp/counts_compiled_hg.tsv
python scripts/rarefy_counts_matrix.py --input tmp/counts_compiled_hg.tsv --n_reads example_data/total_read_count.tsv --subset auxiliary_files/included_args.txt --output rarefied_counts_hg.tsv
rm tmp/counts_compiled_hg.tsv

# Compile and rarefy counts matrix (wastewater samples)
find example_data/mapped_reads_wastewater/filtered -type f > tmp/paths.txt
cat tmp/paths.txt | rev |cut -d '/' -f 1 | rev | cut -d '_' -f 1 > tmp/ids.txt
paste tmp/ids.txt tmp/paths.txt > tmp/sample_paths_ww.txt
rm tmp/paths.txt tmp/ids.txt

python scripts/compile_counts_matrix.py --files tmp/sample_paths_ww.txt --ids auxiliary_files/arg_ids.txt --output tmp/counts_compiled_ww.tsv
python scripts/rarefy_counts_matrix.py --input tmp/counts_compiled_ww.tsv --n_reads example_data/total_read_count.tsv --subset auxiliary_files/included_args.txt --output rarefied_counts_ww.tsv
rm tmp/counts_compiled_ww.tsv

rm -r tmp

# Compile ARG metadata
python scripts/compile_arg_metadata.py --input rarefied_counts_hg.tsv --pathogens auxiliary_files/pathogen_list.txt --arg_index auxiliary_files/arg_index.txt --blastout auxiliary_files/blastout_established.tsv --cluster_dir example_data/clusters --taxonomy_ncbi auxiliary_files/genome_full_lineage.tsv --taxonomy_card auxiliary_files/card_taxonomy.txt --output arg_metadata.tsv
python scripts/get_phylum_distribution.py --input rarefied_counts_hg.tsv --cluster_dir example_data/clusters --taxonomy_ncbi auxiliary_files/genome_full_lineage.tsv --taxonomy_card auxiliary_files/card_taxonomy.txt --output arg_phylum_distribution.tsv
python scripts/get_pathogen_distribution.py --input rarefied_counts_hg.tsv --pathogens auxiliary_files/pathogen_list.txt --cluster_dir example_data/clusters --taxonomy_ncbi auxiliary_files/genome_full_lineage.tsv --taxonomy_card auxiliary_files/card_taxonomy.txt --output arg_pathogen_distribution.tsv