#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --time 8:00:00
#SBATCH -J beta
#SBATCH --mem 64g

source activate qiime2-2022.2

python process.py \
    --tree ${tree} \
    --table-16s ${table_16s} \
    --table-wgs ${table_wgs} \
    --metadata ${metadata} \
    --taxonomy ${taxonomy} \
    --output ${label} \
    --threads 4 \
    --depth-16S ${depth_16S} \
    --depth-WGS ${depth_WGS}
