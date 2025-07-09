### Introduction
This repository contains scripts used to generate the results presented in the paper "Community-promoted antibiotic resistance genes show increased dissemination among pathogens" by Lund el al. 2025. Alongside the scripts, this repository also provides example data which can be used to run the analysis, all files necessary to recreate the main and supplementary figures from the manuscript, as well as a .yml file which can be used to create a conda environment with the required software installed.

### Dependencies
To run the scripts in this repository, the following software are required:
- Python >= 3.12.3
    - biopython >= 1.83
    - numpy >= 1.26.4
    - pandas >= 2.2.2
    - tqdm >= 4.66.4
- Snakemake >= 8.14.0
- R >= 4.3.0
    - data.table >= 1.15.4
    - GGally >= 2.2.1
    - ggpp >= 3.5.1
    - igraph >= 2.0.3
    - patchwork >= 1.2.0
    - scales >= 1.3.0
    - tidyverse >= 2.0.0
- kraken2 >= 2.1.3

A .yml file which can be used to create a conda envioronment with the necessary software installed can be found in /envs.

### Tutorial

Below is a step-by-step guide on how to run the scripts in this repository. Please note that this pipeline is designed to run on Unix based serves.

1. Clone the repository using
    ```
    git clone https://github.com/davidgllund/promoted_args_hg_ww.git
    ```

2. Setup conda environment
    ```
    cd promoted_args_hg_ww
    conda env create -f envs/arg_analysis.yml
    ```

3. To run the analysis pipeline on the provided example data, use the following command:
   ```
    bash scripts/process_diamond_output.sh
    ```

   This will produce the following files:
   - **rarefied_counts_hg.tsv**: Rarefied counts matrix produced from human gut metagenomic samples.
   - **rarefied_counts_ww.tsv**: Rarefied counts matrix produced from wastewater metagenomic samples.
   - **arg_metadata.tsv**: General information about the ARGs included in the analysis.
   - **arg_phylum_distribution.tsv**: Bacterial phyla carrying each of the ARGs included in the analysis.
   - **arg_pathogen_distribution.tsv**: Pathogenis carrying each of the ARGs included in the analysis.
   - **genetic_compatibility.txt**: Euclidean distance between the 5mer distributions of the included ARGs and representative genomes from bacterial pathogens.
  
4. To generate the main and supplementary figures from the paper, use the following command:
   ```
    Rscripts scripts/plots_for_paper.R
    ```

   This will produce a total of 14 files, each named "figure_X.pdf" or "figure_sX.pdf", where X is the corresponding figure number.
