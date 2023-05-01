import click
import bp
import re
import pandas as pd
import pandas.testing as pdt
import skbio
import qiime2


ASV = re.compile('^[ATGC]{90}')


def mark_multifurcations(tree):
    for n in tree.non_tips(include_self=True):
        n.has_multifurcation = 0
    for n in tree.tips():
        if ASV.match(n.name) is not None:
            n.parent.has_multifurcation += 1

def test_mark_multifurcations():
    t = skbio.TreeNode.read(["((A,B)C,(D,E)F)root;"])
    t.find('A').name = 'A' * 100
    t.find('B').name = 'T' * 100
    mark_multifurcations(t)
    assert t.find('C').has_multifurcation == 2
    assert t.find('F').has_multifurcation == 0
test_mark_multifurcations()


def cut(tree):
    mapping = []
    removed = set()
    for n in list(tree.preorder()):
        if not n.is_tip() and n.has_multifurcation > 1:
            name = 'mf-%d' % n.id
            for tip in n.children:
                mapping.append([tip.name, name])
                tip.parent = None
                removed.add(tip)

            n.children = []
            n.name = name
        elif n.is_tip() and n not in removed:
            mapping.append([n.name, n.name])

    return pd.DataFrame(mapping, columns=['feature-id', 'collapsed-tip'])


def test_cut():
    t = skbio.TreeNode.read(["((A,B)C,(D,E)F)root;"])
    t.has_multifurcation = 0
    t.find('F').has_multifurcation = 0
    t.find('C').has_multifurcation = 2
    t.assign_ids()
    exp = pd.DataFrame([['A', 'mf-4'],
                        ['B', 'mf-4'],
                        ['D', 'D'],
                        ['E', 'E']],
                       columns=['feature-id', 'collapsed-tip'])
    obs = cut(t)
    pdt.assert_frame_equal(obs, exp)

    assert t.children[0].name == 'mf-4'
    assert t.children[0].is_tip()
    assert t.children[1].name == 'F'
test_cut()


@click.command()
@click.option('--tree', type=click.Path(exists=True), required=True)
@click.option('--output', type=click.Path(exists=False), required=True)
def process(tree, output):
    tree = bp.to_skbio_treenode(bp.parse_newick(open(tree).read()))

    tree.assign_ids()

    mark_multifurcations(tree)
    mapping = cut(tree)

    # why is this needed?
    for n in tree.traverse(include_self=False):
        if n.length < 0:
            n.length = 0

    mapping.to_csv(output + '.tsv', sep='\t', index=False, header=True)
    tree.write(output + '.nwk')
    ar = qiime2.Artifact.import_data('Phylogeny[Rooted]', tree)
    ar.save(output + '.nwk.qza')

if __name__ == '__main__':
    process()
