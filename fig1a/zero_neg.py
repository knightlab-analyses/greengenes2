import sys
import bp
import qiime2
import numpy as np

t = bp.parse_newick(open(sys.argv[1]).read())
lengths = np.zeros(len(t.B))
for i in range(len(t.B)):
    lengths[i] = t.length(i)
lengths[lengths < 0] = 0.
t.set_lengths(lengths)
sk = bp.to_skbio_treenode(t)
ar = qiime2.Artifact.import_data('Phylogeny[Rooted]', sk)
ar.save(sys.argv[2])

