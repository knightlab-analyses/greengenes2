#!/bin/bash

# note this assumes blast results are already computed and available
sbatch --export phylogeny=../fig1a/2022.10.phylogeny.collapsed_multifurcation.nwk.nwk,collapse_map=../fig1a/2022.10.phylogeny.collapsed_multifurcation.nwk.tsv,blast_results=gg_13_8_aligned.tsv,output=gg_13_8.tsv,label=GG_13_8 submit-collapsed.sbatch
sbatch --export phylogeny=../fig1a/2022.10.phylogeny.collapsed_multifurcation.nwk.nwk,collapse_map=../fig1a/2022.10.phylogeny.collapsed_multifurcation.nwk.tsv,blast_results=silva_138_aligned.tsv,output=silva_138.tsv,label=SILVA_138 submit-collapsed.sbatch
sbatch --export phylogeny=../fig1a/2022.10.phylogeny.collapsed_multifurcation.nwk.nwk,collapse_map=../fig1a/2022.10.phylogeny.collapsed_multifurcation.nwk.tsv,blast_results=gtdb_207_aligned.tsv,output=gtdb_207.tsv,label=GTDB_207 submit-collapsed.sbatch
