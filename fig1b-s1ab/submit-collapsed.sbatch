#!/bin/bash
#SBATCH -J fig1b
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem 24g
#SBATCH --time 24:00:00

source activate qiime2-2022.2

python process-collapsed.py \
    --collapse-map ${collapse_map} \
    --blast-results ${blast_results} \
    --output ${output} \
    --label ${label}

qiime empress tree-plot \
    --i-tree ${phylogeny}.qza \
    --m-feature-metadata-file ${output} \
    --p-no-shear-to-feature-metadata \
    --o-visualization ${output}.qzv
