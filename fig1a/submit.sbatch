#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem 64g
#SBATCH --time 4:00:00
#SBATCH -J cut

source activate qiime2-2022.2

version=2022.10
label=fig1a_multicollapse.tsv

curl --output fig1a.tsv.gz "https://zenodo.org/record/7884635/files/fig1a.tsv.gz?download=1"

python collapse-multifurcation.py \
    --tree ${version}.phylogeny.asv.nwk \
    --output ${version}.phylogeny.collapsed_multifurcation.nwk

python zero_neg.py ${version}.phylogeny.collapsed_multifurcation.nwk.nwk tree.qza
python empress-metadata-collapsed-tree.py \
    --collapse-map ${version}.phylogeny.collapsed_multifurcation.nwk.tsv \
    --tip-association fig1a.tsv.gz \
    --output ${label}

qiime empress tree-plot \
    --i-tree tree.qza \
    --m-feature-metadata-file ${label} \
    --p-no-shear-to-feature-metadata \
    --o-visualization ${label}.qzv
