import skbio
import sys
import math

f = sys.argv[1]
splits = int(sys.argv[2])

data = [(r.metadata['id'], str(r)) for r in skbio.read(f, format='fasta')]
size = math.ceil(len(data) / splits)

for i, chunk in enumerate(range(0, len(data), size)):
    with open(f + '.%d' % i, 'w') as fp:
        for id, seq in data[chunk:chunk + size]:
            fp.write('>%s\n%s\n' % (id, seq))
