#!/bin/bash

curl -O http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza
curl -O http://ftp.microbio.me/greengenes_release/gg_13_8_otus/rep_set/99_otus.fasta
curl -O https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
curl -O https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_all/ssu_all_r207.tar.gz

gunzip SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
tar xzf ssu_all_r207.tar.gz

mkdir -p $PANFS/fast_blast_silva
mkdir -p $PANFS/fast_blast_gtdb
mkdir -p $PANFS/fast_blast_gg138

# the gg2 backbone
unzip 2022.10.backbone.full-length.fna.qza
cp a53d9300-5c5c-4774-a2e8-a5e23904f1ae/data/dna-sequences.fasta* $PANFS/

cp SILVA_138.1_SSURef_NR99_tax_silva.fasta $PANFS/fast_blast_silva/
cp ssu_all_r207.fna $PANFS/fast_blast_gtdb/
cp 99_otus.fasta $PANFS/fast_blast_gg138/

a=$(sbatch --parsable --export f=$PANFS/fast_blast_silva/SILVA_138.1_SSURef_NR99_tax_silva.fasta split_fasta.sbatch)
b=$(sbatch --parsable --export f=$PANFS/fast_blast_gtdb/ssu_all_r207.fna split_fasta.sbatch)
c=$(sbatch --parsable --export f=$PANFS/fast_blast_gg138/99_otus.fasta split_fasta.sbatch)

d=$(sbatch --parsable --dependency=afterok:${a} --export QUERY=$PANFS/fast_blast_silva/SILVA_138.1_SSURef_NR99_tax_silva.fasta,OUTPUT=$PANFS/fast_blast_silva/silva_138_aligned.tsv fast_blast.sbatch)
e=$(sbatch --parsable --dependency=afterok:${b} --export QUERY=$PANFS/fast_blast_gtdb/ssu_all_r207.fna,OUTPUT=$PANFS/fast_blast_gtdb/gtdb_207_aligned.tsv fast_blast.sbatch)
f=$(sbatch --parsable --dependency=afterok:${c} --export QUERY=$PANFS/fast_blast_gg138/99_otus.fasta,OUTPUT=$PANFS/fast_blast_gg138/gg_13_8_aligned.tsv fast_blast.sbatch)

sbatch --dependency=afterok:${d} --export f=$PANFS/fast_blast_silva/silva_138_aligned.tsv merge.sbatch
sbatch --dependency=afterok:${e} --export f=$PANFS/fast_blast_gtdb/gtdb_207_aligned.tsv merge.sbatch
sbatch --dependency=afterok:${f} --export f=$PANFS/fast_blast_gg138/gg_13_8_aligned.tsv merge.sbatch
