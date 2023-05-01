#!/bin/bash

curl -O http://ftp.microbio.me/greengenes_release/current/2022.10.phylogeny.asv.nwk

sbatch submit.sbatch
