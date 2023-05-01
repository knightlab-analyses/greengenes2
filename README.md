Greengenes2 figure generation
-----------------------------

The various figure panels used in the Greengenes2 paper can be constructed from
this repository. 

A few assumptions are made. First, submission scripts are designed for a SLURM
based compute cluser. Second, a QIIME 2 2022.11 environment is available. 
And third, some of the figures depend on `empress` and `iow` which are pip
installable. 

Each directory has a `README.md` and a `submit.sh`. Precomputed data QIIME 2
QZAs are provided for convenience. 

Please note that while we discuss the FINRISK data in the manuscript, we cannot
release any of those data in this repository due to access control 
restrictions.
