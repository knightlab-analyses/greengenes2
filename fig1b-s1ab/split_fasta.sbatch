#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time 24:00:00
#SBATCH -J split
#SBATCH --mem 8g

source activate qiime2-2022.2

python split_fasta.py $f 500
