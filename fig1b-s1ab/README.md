Database comparision
--------------------

Compute for these panels comes in two stages. First, BLAST results must be 
computed and second, the results are summarized.

Processing assumes a high performance storage space is available under $PANFS,
however that location can be readily modified to fit local environments.

First, use `submit_fast_blast.sh` to compute BLAST results. Once computed, then
use `submit-collapsed.sh`. 

The assessment of the BLAST output assumes the multifurcation collapsed map
from figure 1A is available. It also assumes the the following library is
available

```
$ pip install iow
```
