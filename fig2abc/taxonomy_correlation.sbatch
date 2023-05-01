#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time 1:00:00
#SBATCH --mem 8gb

source activate qiime2-2022.2

table_a_label=$(echo ${table_a_label} | tr "_" " ")
table_b_label=$(echo ${table_b_label} | tr "_" " ")
python compare.py \
    --table-a ${table_a} \
    --table-b ${table_b} \
    --output ${output}.extract_label \
    --table-a-name "${table_a_label}" \
    --table-b-name "${table_b_label}" \
    --title ${title} \
    --extract-label

python compare.py \
    --table-a ${table_a} \
    --table-b ${table_b} \
    --output ${output}.no_extract_label \
    --table-a-name "${table_a_label}" \
    --table-b-name "${table_b_label}" \
    --title ${title} 
