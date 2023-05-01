#!/bin/bash

RESULTS=taxonomy-correlation
THDMI_16S=THDMI-16S
THDMI_WGS=THDMI-WGS

mkdir -p $RESULTS

# join strings in an array, see 
# https://stackoverflow.com/a/17841619
function join_by { local IFS="$1"; shift; echo "$*"; }
declare -a jobs
for r in {species,genus,family,order,class,phylum}
do
    table_a=${THDMI_16S}/gg2-feature-table-${r}.qza
    table_a_label="GG2_16S_log10(rel._abund)"
    table_b=${THDMI_WGS}/gg2-feature-table-${r}.qza
    table_b_label="GG2_WGS_log10(rel._abund)"
    output=${RESULTS}/gg2_16s_gg2_wgs_${r}
    jobs+=($(sbatch --export title=${r},table_a=${table_a},table_b=${table_b},output=${output},table_a_label=${table_a_label},table_b_label=${table_b_label} \
        --parsable \
        taxonomy_correlation.sbatch))
done

for r in {species,genus,family,order,class,phylum}
do
    table_a=${THDMI_16S}/gg-feature-table-${r}.qza
    table_a_label="GG138_16S_log10(rel._abund)"
    table_b=${THDMI_WGS}/gg2-feature-table-${r}.qza
    table_b_label="GG2_WGS_log10(rel._abund)"
    output=${RESULTS}/gg138_16s_gg2_wgs_${r}
    jobs+=($(sbatch --export title=${r},table_a=${table_a},table_b=${table_b},output=${output},table_a_label=${table_a_label},table_b_label=${table_b_label} \
        --parsable \
        taxonomy_correlation.sbatch))
done

for r in {genus,family,order,class,phylum}
do
    table_a=${THDMI_16S}/silva-feature-table-${r}.qza
    table_a_label="SILVA_16S_log10(rel._abund)"
    table_b=${THDMI_WGS}/gg2-feature-table-${r}.qza
    table_b_label="GG2_WGS_log10(rel._abund)"
    output=${RESULTS}/silva_16s_gg2_wgs_${r}
    jobs+=($(sbatch --export title=${r},table_a=${table_a},table_b=${table_b},output=${output},table_a_label=${table_a_label},table_b_label=${table_b_label} \
        --parsable \
        taxonomy_correlation.sbatch))
done

for r in {species,genus,family,order,class,phylum}
do
    table_a=${THDMI_16S}/gg2nb-feature-table-${r}.qza
    table_a_label="GG2nb_16S_log10(rel._abund)"
    table_b=${THDMI_WGS}/gg2-feature-table-${r}.qza
    table_b_label="GG2_WGS_log10(rel._abund)"
    output=${RESULTS}/gg2nb_16s_gg2_wgs_${r}
    jobs+=($(sbatch --export title=${r},table_a=${table_a},table_b=${table_b},output=${output},table_a_label=${table_a_label},table_b_label=${table_b_label} \
        --parsable \
        taxonomy_correlation.sbatch))
done

dependency=$(join_by : ${jobs[@]})
sbatch --dependency=afterok:${dependency} \
    --export results=$RESULTS \
    per-sample.sbatch
