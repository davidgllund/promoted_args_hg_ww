#!/usr/bin/bash 
if [[ ! -d bacterial_taxonomy ]]; then
    mkdir bacterial_taxonomy
fi

while read line; do

    kraken2 --paired "$line"_1_QC.fastq.gz "$line"_2_QC.fastq.gz --out kraken_output.csv --db /storage/dlund/Michaela_reloaded/data/kraken/kraken_db/ --threads 65 --gzip-compressed --confidence 0.1 --use-names

    python /home/dlund/Michaela_reloaded/scripts/count_bacteria.py kraken_output.csv

    rm kraken_output.csv
done < $1
echo "Finished all"

echo "Start merge in python"
python3 /home/dlund/Michaela_reloaded/scripts/merge_taxonomy.py
wait

echo "DONE!"
