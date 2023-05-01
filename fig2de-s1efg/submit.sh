#!/bin/bash

dm1=thdmi-gg2-table_16S_even_9000.weighted_unifrac.qza
dm2=thdmi-gg2-table_WGS_even_2000000.weighted_unifrac.qza
dm1label="GG2-16S-Weighted-UniFrac"
dm2label="GG2-WGS-Weighted-UniFrac"
metadata=thdmi-metadata.tsv
output=fig2e
sbatch --export dm1=${dm1},dm2=${dm2},dm1label=${dm1label},dm2label=${dm2label},metadata=${metadata},output=${output} submit-effect-size.sbatch

dm1=thdmi-silva-table_16S_even_9000.weighted_unifrac.qza
dm2=thdmi-gg2-table_WGS_even_2000000.weighted_unifrac.qza
dm1label="SILVA-16S-Weighted-UniFrac"
dm2label="GG2-WGS-Weighted-UniFrac"
metadata=thdmi-metadata.tsv
output=figs1e
sbatch --export dm1=${dm1},dm2=${dm2},dm1label=${dm1label},dm2label=${dm2label},metadata=${metadata},output=${output} submit-effect-size.sbatch

dm1=thdmi-na-table_16S_even_9000.braycurtis.qza
dm2=thdmi-na-table_WGS_even_2000000.braycurtis.qza
dm1label="GG2-16S-Bray-Curtis"
dm2label="GG2-WGS-Bray-Curtis"
metadata=thdmi-metadata.tsv
output=fig2d
sbatch --export dm1=${dm1},dm2=${dm2},dm1label=${dm1label},dm2label=${dm2label},metadata=${metadata},output=${output} submit-effect-size.sbatch
