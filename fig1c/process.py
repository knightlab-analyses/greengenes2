import pandas as pd
import matplotlib.pyplot as plt
import bp
import click
import biom
import seaborn as sns
import numpy as np
import re
import skbio


ASV = re.compile('^[ATGC]{90}')

# test in ../collapse-multifurcation/process.py
def mark_multifurcations(tree):
    for n in tree.non_tips(include_self=True):
        n.has_multifurcation = 0
    for n in tree.tips():
        if ASV.match(n.name) is not None:
            n.parent.has_multifurcation += 1


def get_branch_length(tree, table):
    in_table = set(table.ids(axis='observation'))
    bl = {}

    tree_total = 0.
    for n in tree.traverse(include_self=False):
        if n.is_tip() and ASV.match(n.name) is not None:
            continue
        else:
            if not n.is_tip() and n.has_multifurcation > 1:
                continue
            tree_total += n.length

    for tip in tree.tips():
        if tip.name in in_table:
            bl[tip.name] = tip.length
            if tip.parent.has_multifurcation > 1:
                bl[tip.name] += tip.parent.length

    table = table.filter(set(bl), axis='observation',
                         inplace=False).remove_empty()
    table_ids = table.ids(axis='observation')

    result = {}
    for v, i, m in table.iter(dense=False):
        total_length = 0.
        for asv_idx in v.indices:
            total_length += bl[table_ids[asv_idx]]
        result[i] = total_length / tree_total

    return result


def test_get_branch_length():
    tree = skbio.TreeNode.read(['(((a:0.1,b:0.2)c:0.1,((x:0.3,y:0.4)z:0.1)w:0.5)cool:0.1);'])
    mark_multifurcations(tree)
    table = biom.Table(np.array([[0, 1, 2],
                                 [3, 0, 4]]),
                       ['x', 'y'],
                       ['S1', 'S2', 'S3'])
    exp = {'S1': 0.4 / 1.8,
           'S2': 0.3 / 1.8,
           'S3': 0.7 / 1.8}
    obs = get_branch_length(tree, table)
    assert obs.keys() == exp.keys()
    assert np.allclose(list(obs.values()), list(exp.values()))

test_get_branch_length()


@click.command()
@click.option('--tree', type=click.Path(exists=True), required=True)
@click.option('--table', type=click.Path(exists=True), required=True)
@click.option('--metadata', type=click.Path(exists=True), required=True)
@click.option('--output', type=click.Path(exists=False), required=True)
def summarize(tree, table, metadata, output):
    tree = bp.parse_newick(open(tree).read())
    tips = {tree.name(i) for i in range(len(tree.B) - 1)
            if tree.B[i] and not tree.B[i+1]}
    md = pd.read_csv(metadata, sep='\t').set_index('#SampleID')

    counts = md['empo_3'].value_counts()
    keep = counts[counts > 50]
    md = md[md['empo_3'].isin(keep.index)]

    table = biom.load_table(table)
    overlap = set(table.ids(axis='observation')) & tips
    table = table.filter(overlap, axis='observation')
    table = table.subsample(1000)

    tree = tree.shear(overlap).collapse()
    tree = bp.to_skbio_treenode(tree)
    mark_multifurcations(tree)

    def transform(v, i, m):
        vnorm = v / v.sum()
        v[vnorm < 0.01] = 0
        return v
    table = table.transform(transform).remove_empty()

    overlap = set(table.ids()) & set(md.index)
    table = table.filter(overlap).remove_empty()
    md = md.loc[list(overlap)]

    bl = get_branch_length(tree, table)

    # we may loose samples with low depth and features not in the tree
    md = md.loc[list(bl.keys())]

    results = []
    for env, env_grp in md.groupby('empo_3'):
        for sample in env_grp.itertuples():
            v = bl[sample.Index]
            if v > 0:
                # v = np.log10(v)
                results.append((sample.Index, env, v))

    results = pd.DataFrame(results, columns=['#SampleID', 'EMPO3', 'Branch Length'])
    results.to_csv(output + '.tsv', sep='\t', index=False, header=True)

    plt.figure(figsize=(12,8))
    sns.set_style('white')
    sns.set(font_scale=2)
    medians = results.groupby('EMPO3')['Branch Length'].median()
    medians = medians.sort_values()
    #sns.barplot(data=results, x='EMPO3', y='log10(Branch Length)', order=list(medians.index))
    sns.swarmplot(data=results, x='EMPO3', y='Branch Length', order=list(medians.index), size=1)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output + '.pdf')


if __name__ == '__main__':
    summarize()
