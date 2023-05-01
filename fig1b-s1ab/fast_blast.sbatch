#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem 8gb
#SBATCH --time 12:00:00
#SBATCH -J blasting
#SBATCH -a 1-500

source activate qiime2-2022.2

if [[ ! -f ${QUERY}.${SLURM_ARRAY_TASK_ID} ]];
then
    exit 0
fi

blastn \
    -query ${QUERY}.${SLURM_ARRAY_TASK_ID} \
    -db $PANFS/dna-sequences.fasta \
    -out ${OUTPUT}.${SLURM_ARRAY_TASK_ID} \
    -outfmt 6 \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -num_threads 1
