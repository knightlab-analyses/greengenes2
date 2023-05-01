#!/bin/bash
#SBATCH -J fig1c
#SBATCH -N 1 
#SBATCH -c 1
#SBATCH --mem 20g
#SBATCH --time 24:00:00

source activate qiime2-2022.2

tree=2022.10.phylogeny.asv.nwk

ctx_90=Deblur_2021.09-Illumina-16S-V4-90nt-dd6875

cut -f 1 emp_table_90.tsv | grep -v "^#" | \
    redbiom fetch samples \
        --context ${ctx_90} \
        --output emp_table_90.biom 

python process.py \
    --tree ${tree} \
    --metadata emp_table_90.tsv \
    --table emp_table_90.biom \
    --output emp90
