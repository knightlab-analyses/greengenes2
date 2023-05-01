import click
import pandas as pd
import bp
import io


def parse_blast(f):
    obs = set()
    for line in f:
        parts = line.split('\t')
        subject = parts[1]
        obs.add(subject)
    return obs


def test_parse_blast():
    data = """1111886	AB258386	97.467	1382	33	2	1	1381	10	1390	0.0	2357
1111883	G003221035	96.823	1196	30	5	1	1188	28	1223	0.0	1991
1111882	MJ034-2-barcode65-umi2612bins-ubs-156	99.041	1460	13	1	21	1479	1	1460	0.0	2617
1111879	G000016665"""
    exp = {'AB258386', 'G003221035', 'MJ034-2-barcode65-umi2612bins-ubs-156',
           'G000016665'}
    obs = parse_blast(io.StringIO(data))
    assert obs == exp
test_parse_blast()


@click.command()
@click.option('--collapse-map', type=click.Path(exists=True), required=True)
@click.option('--blast-results', type=click.Path(exists=True), required=True)
@click.option('--output', type=click.Path(exists=False), required=True)
@click.option('--label', type=str, required=True)
def fig1bcd(collapse_map, blast_results, output, label):
    collapse_map = pd.read_csv(collapse_map, sep='\t').set_index('feature-id')

    observed = parse_blast(open(blast_results))
    feature_md = []

    for k, kgrp in collapse_map.groupby('collapsed-tip'):
        feature_md.append([k, any([v in observed for v in kgrp.index])])

    feature_md = pd.DataFrame(feature_md, columns=['FeatureID', label])
    feature_md.to_csv(output, sep='\t', index=False, header=True)


if __name__ == '__main__':
    fig1bcd()
