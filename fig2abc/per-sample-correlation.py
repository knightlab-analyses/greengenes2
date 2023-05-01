import re
import pandas as pd
import glob
import click
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
import pandas.testing as pdt
import os


# gg2_16s_wol2_wgs_species.no_extract_label.tsv
EXTRACT = re.compile(r"^(\w+_\w+)_(\w+_\w+)_(\w+)\.(\w+)\.tsv$")
def extract(filename):
    match = EXTRACT.match(filename)
    if match is None:
        raise ValueError("cannot parse %s" % filename)
    else:
        method_a, method_b, rank, extract = match.groups()
        return (method_a, method_b, rank, extract)


def test_extract():
    tests = [('gg2_16s_wol2_wgs_species.no_extract_label.tsv',
               ('gg2_16s', 'wol2_wgs', 'species', 'no_extract_label')),]
    for test, exp in tests:
        obs = extract(test)
        assert obs == exp
test_extract()


def concat_dataframes(directory):
    extra_cols = {'sample-id', 'feature', 'method-a', 'method-b', 'rank', 'extract'}
    dfs = []
    for f in glob.glob(directory + '/*.tsv'):
        fname = os.path.basename(f)
        if fname.startswith('mannu') or fname.startswith('ttest'):
            continue
        method_a, method_b, rank, extracted_names = extract(fname)
        df = pd.read_csv(f, sep='\t')
        labels = [c for c in df.columns if c not in extra_cols]
        assert len(labels) == 2
        method_name_a, method_name_b = labels
        df.rename(columns={method_name_a: 'method_a',
                           method_name_b: 'method_b'}, inplace=True)
        df['method-a'] = method_name_a
        df['method-b'] = method_name_b
        df['rank'] = rank
        df['extract'] = extracted_names
        dfs.append(df)
    return pd.concat(dfs)


def compress(df):
    extra_cols = ['method-a', 'method-b', 'rank', 'extract']

    # there should only be two columns remaining, if not we have a problem
    rel_a, rel_b = (set(df.columns) - set(extra_cols)) - {'sample-id', 'feature'}
    per_sample = []
    for (sample, ma, mb, rank, ext), grp in df.groupby(['sample-id'] + extra_cols):
        if len(grp) < 2:
            # need at least 2 entries for a correlation
            continue
        r, _ = ss.pearsonr(grp[rel_a], grp[rel_b])
        p, _ = ss.spearmanr(grp[rel_a], grp[rel_b])
        per_sample.append((sample, r, p, ma, mb, rank, ext))
    return pd.DataFrame(per_sample,
                        columns=['sample-id',
                                 'Pearson correlation',
                                 'Spearman correlation'] + extra_cols)


def test_compress():
    df = pd.DataFrame([['x', 'foo', 0.1, 0.1, 'a', 'b', 'c', 'd'],
                       ['x', 'bar', 0.5, 0.2, 'a', 'b', 'c', 'd'],
                       ['x', 'baz', 0.4, 0.7, 'a', 'b', 'c', 'd'],
                       ['y', 'bar', 0.6, 0.3, 'a', 'b', 'c', 'd'],
                       ['y', 'baz', 0.4, 0.7, 'a', 'b', 'c', 'd']],
                      columns=['sample-id', 'feature', 'X', 'Y', 'method-a',
                               'method-b', 'rank', 'extract'])
    exp = pd.DataFrame([['x',
                         ss.pearsonr([0.1, 0.5, 0.4], [0.1, 0.2, 0.7])[0],
                         ss.spearmanr([0.1, 0.5, 0.4], [0.1, 0.2, 0.7])[0],
                         'a', 'b', 'c', 'd'],
                        ['y',
                         ss.pearsonr([0.6, 0.4], [0.3, 0.7])[0],
                         ss.spearmanr([0.6, 0.4], [0.3, 0.7])[0],
                         'a', 'b', 'c', 'd']],
                       columns=['sample-id', 'Pearson correlation',
                                'Spearman correlation', 'method-a',
                                'method-b', 'rank', 'extract'])
    obs = compress(df)
    pdt.assert_frame_equal(obs, exp)
test_compress()

@click.command()
@click.option('--directory', type=click.Path(exists=True), required=True)
@click.option('--output-directory', type=click.Path(exists=True), required=True)
def per_sample_correlation(directory, output_directory):
    df = concat_dataframes(directory)
    df_compressed = compress(df)

    order = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    df_compressed['method-a'] = df_compressed['method-a'].apply(lambda x: x.split(' log')[0])
    df_compressed['method-b'] = df_compressed['method-b'].apply(lambda x: x.split(' log')[0])
    sns.set(font_scale=2)
    sns.set_style(style='white')

    nonblank = {i for i in df_compressed['sample-id'] if not 'blank' in i.lower()}
    df_compressed = df_compressed[df_compressed['sample-id'].isin(nonblank)]

    remap = {'GG2 16S': 'GG2 Phylogenetic taxonomy',
             'GG2nb 16S': 'GG2 Naive Bayes',
             'GTDB': 'GTDB r207',
             'SILVA 16S': 'SILVA 138 Naive Bayes',
             'GG138 16S': 'GG 13_8 Naive Bayes'}
    df_compressed['method-a'] = df_compressed['method-a'].apply(lambda x: remap.get(x, x))
    df_compressed['method-b'] = df_compressed['method-b'].apply(lambda x: remap.get(x, x))

    comparisons = [
        ('gg2nb_16s_v_gg2_wgs', ('GG2 Naive Bayes', 'GG2 WGS')),
        ('gg2_16s_v_gg2_wgs', ('GG2 Phylogenetic taxonomy', 'GG2 WGS')),
        ('silva_16s_v_gg2_wgs', ('SILVA 138 Naive Bayes', 'GG2 WGS')),
        ('gg_16s_v_gg2_wgs', ('GG 13_8 Naive Bayes', 'GG2 WGS'))
    ]

    for name, (ma, mb) in comparisons:
        subset = df_compressed[(df_compressed['method-a'] == ma) &
                               (df_compressed['method-b'] == mb)]

        for extract, grp in subset.groupby('extract'):
            print(name, extract)
            plt.figure(figsize=(16, 8))
            sns.boxplot(data=grp, x='rank', y='Pearson correlation', hue='method-a', order=order)
            ax = plt.gca()
            ax.set_xlabel('')
            ax.set_title('Taxonomy Correlation', fontsize=24)
            plt.savefig(output_directory + '/pearson_%s_%s.pdf' % (name, extract))
            plt.figure(figsize=(16, 8))

            plt.figure(figsize=(8, 6))
            percs_25 = []
            percs_50 = []
            percs_75 = []

            if 'SILVA' in ma or 'SILVA' in mb:
                ranks = ['class', 'order', 'family', 'genus']
                xpos = [1, 2, 3, 4]
                xlabels = ['Class', 'Order', 'Family', 'Genus']
            else:
                ranks = ['class', 'order', 'family', 'genus', 'species']
                xpos = [1, 2, 3, 4, 5]
                xlabels = ['Class', 'Order', 'Family', 'Genus', 'Species']

            for label in ranks:
                subset = grp[grp['rank'] == label]['Pearson correlation'].describe()
                percs_25.append(subset['25%'])
                percs_50.append(subset['50%'])
                percs_75.append(subset['75%'])
                print("%s vs %s; %s: %0.2f" % (ma, mb, label, subset['50%']))
            plt.plot(xpos, percs_50, color='b')
            plt.plot(xpos, percs_25, linestyle='--', color='b')
            plt.plot(xpos, percs_75, linestyle='--', color='b')
            ax = plt.gca()
            ax.set_ylim(0, 1)
            ax.set_xlim(1, max(xpos))
            ax.set_xticks(xpos)
            ax.set_xticklabels(xlabels)
            ax.set_xlabel('Pearson correlation')
            ax.set_title("%s vs. %s" % (ma, mb))
            plt.tight_layout()
            plt.savefig(output_directory + '/pearson_lineplot_%s_%s.pdf' % (name, extract))

    subset = df_compressed[(df_compressed['method-a'] == 'GG2 Phylogenetic taxonomy') &
                           (df_compressed['method-b'] == 'GG2 WGS') &
                           (df_compressed['extract'] == 'no_extract_label')]
    print(subset.shape)
    plt.figure(figsize=(8, 6))
    percs_25 = []
    percs_50 = []
    percs_75 = []
    for label in ['class', 'order', 'family', 'genus', 'species']:
        subset2 = subset[subset['rank'] == label]['Pearson correlation'].describe()
        percs_25.append(subset2['25%'])
        percs_50.append(subset2['50%'])
        percs_75.append(subset2['75%'])
    plt.plot([1, 2, 3, 4, 5], percs_50, color='b')
    plt.plot([1, 2, 3, 4, 5], percs_25, linestyle='--', color='b')
    plt.plot([1, 2, 3, 4, 5], percs_75, linestyle='--', color='b')

    subset = df_compressed[(df_compressed['method-a'] == 'GG2 Naive Bayes') &
                           (df_compressed['method-b'] == 'GG2 WGS') &
                           (df_compressed['extract'] == 'no_extract_label')]
    print(subset.shape)
    percs_25 = []
    percs_50 = []
    percs_75 = []
    for label in ['class', 'order', 'family', 'genus', 'species']:
        subset2 = subset[subset['rank'] == label]['Pearson correlation'].describe()
        percs_25.append(subset2['25%'])
        percs_50.append(subset2['50%'])
        percs_75.append(subset2['75%'])
    plt.plot([1, 2, 3, 4, 5], percs_50, color='r')
    plt.plot([1, 2, 3, 4, 5], percs_25, linestyle='--', color='r')
    plt.plot([1, 2, 3, 4, 5], percs_75, linestyle='--', color='r')

    ax = plt.gca()
    ax.set_ylim(0, 1)
    ax.set_xlim(1, 5)
    ax.set_xticklabels(['Class', 'Order', 'Family', 'Genus', 'Species'])
    ax.set_xlabel('Pearson correlation')
    ax.set_title('GG2 16S vs. GG2 WGS')
    plt.tight_layout()
    plt.savefig(output_directory + '/pearson_lineplot_combo.pdf')


if __name__ == '__main__':
    per_sample_correlation()
