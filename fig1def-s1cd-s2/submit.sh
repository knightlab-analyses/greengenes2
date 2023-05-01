#!/bin/bash

curl -O http://ftp.microbio.me/greengenes_release/2022.10/2022.10.taxonomy.asv.tsv.qza
curl -O http://ftp.microbio.me/greengenes_release/2022.10/2022.10.phylogeny.asv.nwk.qza

version=2022.10
taxonomy=${version}.taxonomy.asv.tsv.qza
tree=${version}.phylogeny.asv.nwk.qza

depth_16S=10000
depth_WGS=1000000
metadata=thdmi-metadata.txt
table_16s=thdmi-sepp-feature-table.biom.qza
table_wgs=thdmi-wol2-feature-table.biom.qza
label=THDMI
sbatch --export depth_16S=${depth_16S},depth_WGS=${depth_WGS},tree=${tree},taxonomy=${taxonomy},metadata=${metadata},table_16s=${table_16s},table_wgs=${table_wgs},label=${label} submit.sbatch

depth_16S=1000
depth_WGS=50000
metadata=emp500-metadata.txt
table_16s=emp500-sepp-feature-table.biom.qza
table_wgs=emp500-wol2-feature-table.biom.qza
label=EMP500
sbatch --export depth_16S=${depth_16S},depth_WGS=${depth_WGS},tree=${tree},taxonomy=${taxonomy},metadata=${metadata},table_16s=${table_16s},table_wgs=${table_wgs},label=${label} submit.sbatch
