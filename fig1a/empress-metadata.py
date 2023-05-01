import click
from qiime2.plugins import empress
import pandas as pd
import biom

def rewrite_agp_table(t, seqs):
    f = open(seqs)
    id = iter(f)
    seq = iter(f)
    mapping = {}
    for i, s in zip(id, seq):
        i = i[1:-1]
        s = s[:-1]
        mapping[i] = s
    return t.update_ids(mapping, axis='observation')

@click.command()
@click.option('--emp-table-150', type=click.Path(exists=True), required=True)
@click.option('--emp-table-125', type=click.Path(exists=True), required=True)
@click.option('--emp-table-100', type=click.Path(exists=True), required=True)
@click.option('--emp-table-90', type=click.Path(exists=True), required=True)
@click.option('--agp-table', type=click.Path(exists=True), required=True)
@click.option('--agp-seqs', type=click.Path(exists=True), required=True)
@click.option('--gg-taxonomy', type=click.Path(exists=True), required=True)
@click.option('--output', type=click.Path(exists=False), required=True)
def fig1a(emp_table_150, emp_table_125, emp_table_100, emp_table_90, agp_table, agp_seqs, gg_taxonomy, output):
    emp_150 = biom.load_table(emp_table_150)
    emp_125 = biom.load_table(emp_table_125)
    emp_100 = biom.load_table(emp_table_100)
    emp_90 = biom.load_table(emp_table_90)
    agp = biom.load_table(agp_table)
    agp = rewrite_agp_table(agp, agp_seqs)

    agp_ids_150 = set(agp.ids(axis='observation'))
    agp_ids_125 = {i[:125] for i in agp_ids_150}
    agp_ids_100 = {i[:100] for i in agp_ids_150}
    agp_ids_90 = {i[:90] for i in agp_ids_150}
    emp_ids_150 = set(emp_150.ids(axis='observation'))
    emp_ids_125 = set(emp_125.ids(axis='observation'))
    emp_ids_100 = set(emp_100.ids(axis='observation'))
    emp_ids_90 = set(emp_90.ids(axis='observation'))

    tax = pd.read_csv(gg_taxonomy, sep='\t').set_index('Feature ID')
    tax['phylum'] = tax['Taxon'].apply(lambda x: x.split('; ')[1])
    tax['phylum'] = tax['phylum'].apply(lambda x: 'Other' if x == 'p__' else x)
    top20 = {k: k for k in tax['phylum'].value_counts().head(20).index}
    tax['phylum_top20'] = tax['phylum'].apply(lambda x: top20.get(x, 'Other'))

    lengths = (150, 125, 100, 90)
    emp_ids = (emp_ids_150, emp_ids_125, emp_ids_100, emp_ids_90)
    agp_ids = (agp_ids_150, agp_ids_125, agp_ids_100, agp_ids_90)

    feature_md = []
    for r in tax.itertuples():
        is_emp = any([r.Index[:length] in s
                      for length, s in zip(lengths, emp_ids)])
        is_agp = any([r.Index[:length] in s
                      for length, s in zip(lengths, agp_ids)])
        feature_md.append([r.Index, is_emp, is_agp, r.phylum_top20])

    feature_md = pd.DataFrame(feature_md, columns=['feature-id', 'is_emp', 'is_agp', 'phylum_top20'])
    feature_md.to_csv(output, sep='\t', index=False, header=True)


if __name__ == '__main__':
    fig1a()
