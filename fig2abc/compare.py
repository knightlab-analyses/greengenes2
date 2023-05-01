import click
import qiime2
import biom
import numpy as np
import pandas.testing as pdt
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
import pandas as pd
from gg2.species_report import extract_label


def test_overlap_tables():
    a = biom.Table(np.array([[0, 1, 2],
                             [3, 4, 5],
                             [6, 7, 8]]),
                   ['a', 'b', 'c'],
                   ['x', 'y', 'z'])
    b = biom.Table(np.array([[0, 1, 2],
                             [3, 4, 5]]),
                   ['a', 'b'],
                   ['x', 'y', 'w'])
    exp_a = biom.Table(np.array([[0, 1],
                                 [3, 4],
                                 [6, 7]]),
                       ['a', 'b', 'c'],
                       ['x', 'y'])
    exp_b = biom.Table(np.array([[0, 1],
                                 [3, 4]]),
                       ['a', 'b'],
                       ['x', 'y'])

    obs_a, obs_b = overlap_tables(a, b)
    assert obs_a == exp_a
    assert obs_b == exp_b


def overlap_tables(tab_a, tab_b, extract=True):
    def collapse_extract(i, m):
        # account for polyphyletic labelings
        return extract_label(i.split(';')[-1])

    def collapse_noextract(i, m):
        return i.split(';')[-1]

    if extract:
        collapse_f = collapse_extract
    else:
        collapse_f = collapse_noextract

    tab_a = tab_a.collapse(collapse_f, axis='observation', norm=False)
    tab_b = tab_b.collapse(collapse_f, axis='observation', norm=False)
    tab_a.del_metadata()
    tab_b.del_metadata()

    samp_ids = set(tab_a.ids()) & set(tab_b.ids())

    new_a = tab_a.filter(samp_ids, inplace=False).remove_empty()
    new_b = tab_b.filter(samp_ids, inplace=False).remove_empty()

    # regather IDs for sort order
    samp_ids = set(new_a.ids()) & set(new_b.ids())
    samp_ids = sorted(samp_ids)

    # set a minimum for comparison, this makes the lower bound of relative
    # abundance common between the samples
    mins = {i: max(a, b) for i, a, b in zip(new_a.ids(), new_a.min('sample'),
                                            new_b.min('sample'))}
    class minmax:
        def __init__(self):
            self.dropped = 0
        def __call__(self, v, i, m):
            total = v.sum()
            m = mins[i]
            v[v < m] = 0
            self.dropped += (total - v.sum())
            return v

    mm = minmax()
    new_a.transform(mm, inplace=True)
    mm = minmax()
    new_b.transform(mm, inplace=True)

    return (new_a.sort_order(samp_ids),
            new_b.sort_order(samp_ids))
test_overlap_tables()


def test_melt():
    tab_a = biom.Table(np.array([[0, 0, 2],
                                 [3, 4, 5],
                                 [6, 7, 8]]),
                       ['a', 'b', 'c'],
                       ['x', 'y', 'z'])
    tab_b = biom.Table(np.array([[0, 1, 2],
                                 [3, 0, 5],
                                 [9, 7, 8]]),
                       ['a', 'b', 'c'],
                       ['x', 'y', 'z'])
    exp = pd.DataFrame([['x', 'b', 3., 3.],
                        ['x', 'c', 6., 9.],
                        ['y', 'a', 0., 1.],
                        ['y', 'b', 4., 0.],
                        ['y', 'c', 7., 7.],
                        ['z', 'a', 2., 2.],
                        ['z', 'b', 5., 5.],
                        ['z', 'c', 8., 8.]],
                       columns=['sample-id',
                                'feature',
                                'table-A-relative-abundance',
                                'table-B-relative-abundance'])
    obs = melt(tab_a, tab_b)
    obs = obs.sort_values(['sample-id', 'feature'])
    obs.index = list(range(len(obs)))
    pdt.assert_frame_equal(obs, exp)


def melt(tab_a, tab_b):
    feats_a = tab_a.ids(axis='observation')
    feats_b = tab_b.ids(axis='observation')
    observed = []

    # make sure tables are in order
    assert (tab_a.ids() == tab_b.ids()).all()

    for i, v_a, v_b in zip(tab_a.ids(),
                           tab_a.iter_data(dense=False),
                           tab_b.iter_data(dense=False)):

        # construct feature value lookups
        lookup_a = {feats_a[i]: v for i, v in zip(v_a.indices,
                                                  v_a.data)}
        lookup_b = {feats_b[i]: v for i, v in zip(v_b.indices,
                                                  v_b.data)}

        # for each feature, store the relative abundances
        for f in set(lookup_a) | set(lookup_b):
            v1 = lookup_a.get(f, 0.)
            v2 = lookup_b.get(f, 0.)

            observed.append((i,
                             f,
                             v1,
                             v2))

    return pd.DataFrame(observed,
                        columns=['sample-id',
                                 'feature',
                                 'table-A-relative-abundance',
                                 'table-B-relative-abundance'])
test_melt()


def test_norm_if_needed():
    tab_a = biom.Table(np.array([[0, 0, 2],
                                 [3, 4, 5],
                                 [6, 7, 8]]),
                       ['a', 'b', 'c'],
                       ['x', 'y', 'z'])
    tab_b = biom.Table(np.array([[0, 0, .2],
                                 [.3, .4, .5],
                                 [.6, .4, .4]]),
                       ['a', 'b', 'c'],
                       ['x', 'y', 'z'])
    exp_a = biom.Table(np.array([[0, 0, 2 / 15],
                                 [3 / 9, 4 / 11, 5 / 15],
                                 [6 / 9, 7 / 11, 8 / 15]]),
                       ['a', 'b', 'c'],
                       ['x', 'y', 'z'])
    exp_b = tab_b.copy()
    obs_a = norm_if_needed(tab_a)
    obs_b = norm_if_needed(tab_b)
    assert obs_a == exp_a
    assert obs_b == exp_b


def norm_if_needed(tab):
    if (tab.matrix_data.data[:100] > 1).any():
        return tab.norm(inplace=False)
    else:
        return tab


@click.command()
@click.option('--table-a', type=click.Path(exists=True), required=True)
@click.option('--table-b', type=click.Path(exists=True), required=True)
@click.option('--output', type=click.Path(exists=False), required=True)
@click.option('--table-a-name', type=str, required=True)
@click.option('--table-b-name', type=str, required=True)
@click.option('--extract-label', is_flag=True, required=True, default=False)
@click.option('--title', type=str, required=True)
def compare(table_a, table_b, output, table_a_name, table_b_name, extract_label, title):
    table_a = qiime2.Artifact.load(table_a).view(biom.Table)
    table_b = qiime2.Artifact.load(table_b).view(biom.Table)

    table_a = norm_if_needed(table_a)
    table_b = norm_if_needed(table_b)
    table_a, table_b = overlap_tables(table_a, table_b, extract_label)

    df = melt(table_a, table_b)
    df.rename(columns={'table-A-relative-abundance': table_a_name,
                       'table-B-relative-abundance': table_b_name},
              inplace=True)

    n_features = len(df['feature'].unique())

    print("Number of unique features compared: %d" % n_features)
    df.to_csv(output + '.tsv', sep='\t', index=False, header=True)
    plt.figure(figsize=(10,10))
    plt.scatter(df[table_a_name], df[table_b_name], s=1)
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')

    r, p = ss.pearsonr(df[table_a_name], df[table_b_name])

    # per Jamie's ASD paper
    min_ = 1e-4
    max_ = 1

    ax.set_xlim((min_, max_))
    ax.set_ylim((min_, max_))
    ax.text(min_ * 2, 0.6, "r=%0.2f, p=%0.2e" % (r, p), fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.set_xlabel(table_a_name, fontsize=16)
    ax.set_ylabel(table_b_name, fontsize=16)
    plt.tight_layout()
    plt.savefig(output + '.png')


if __name__ == '__main__':
    compare()
