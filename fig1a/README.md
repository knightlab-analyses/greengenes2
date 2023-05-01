Phylogeny big picture
---------------------

This code assumes a standard QIIME 2 environment, with the following additional
packages installed:

```bash
$ pip install empress
$ pip install iow
```

The Empress metadata and the script to produce it are provided. The raw inputs
are not provided to avoid bloat in the repository. Those inputs can be 
reproduced from public data though. 

```bash
$ ctx_150=Deblur_2021.09-Illumina-16S-V4-150nt-ac8c0b
$ ctx_125=Deblur_2021.09-Illumina-16S-V4-125nt-92f954
$ ctx_100=Deblur_2021.09-Illumina-16S-V4-100nt-50b3a2
$ ctx_90=Deblur_2021.09-Illumina-16S-V4-90nt-dd6875
 
$ redbiom search metadata \
>    "where qiita_study_id==10317 and env_package=='human-gut'" | \
>    redbiom fetch samples \
>        --context $ctx_90 \
>        --resolve-ambiguities merge \
>        --output agp_table.biom
$ biom table-ids -i agp_table.biom --observations | \
>    awk '{ print ">" $1 "\n" $1 }' > agp_seqs.fna
$ agp_table=agp_table.biom
$ agp_seqs=agp_seqs.fna

# ftp.microbio.me/emp/release1/otu_tables/deblur/emp_deblur_90bp.release1.biom
$ emp_deblur_90bp.release1.biom

# see http://ftp.microbio.me/greengenes_release/2022.10/
$ gg_taxonomy=2022.10.taxonomy.asv.tsv.gz
$ phylogeny=2022.10.phylogeny.id.nwk

$ biom table-ids -i ${emp_table_90_raw} | \
>     redbiom fetch samples \
>         --context ${ctx_150} \
>         --output emp_table_150.biom &
 
$ biom table-ids -i ${emp_table_90_raw} | \
>     redbiom fetch samples \
>         --context ${ctx_125} \
>         --output emp_table_125.biom &
 
$ biom table-ids -i ${emp_table_90_raw} | \
>     redbiom fetch samples \
>         --context ${ctx_100} \
>         --output emp_table_100.biom &
 
$ biom table-ids -i ${emp_table_90_raw} | \
>     redbiom fetch samples \
>         --context ${ctx_90} \
>         --output emp_table_90.biom &
$ wait

$ python empress-metadata.py \
>     --agp-table ${agp_table} \
>     --agp-seqs ${agp_seqs} \
>     --emp-table-150 emp_table_150.biom \
>     --emp-table-125 emp_table_125.biom \
>     --emp-table-100 emp_table_100.biom \
>     --emp-table-90 emp_table_90.biom \
>     --gg-taxonomy ${gg_taxonomy} \
>     --output fig1a.tsv
```
