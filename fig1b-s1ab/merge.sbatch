#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time 4:00:00
#SBATCH -J merge

cat ${f}.* > ${f}
cp ${f} .
