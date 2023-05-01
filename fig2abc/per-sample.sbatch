#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time 1:00:00
#SBATCH --mem 64gb

source activate qiime2-2022.2

python per-sample-correlation.py \
    --directory ${results} \
    --output-directory ${results}
