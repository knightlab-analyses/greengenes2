import pandas as pd
import qiime2
import click
import seaborn as sns
import scipy.stats as ss
import matplotlib.pyplot as plt


@click.command()
@click.option('--dm1', type=click.Path(exists=True), required=True)
@click.option('--dm2', type=click.Path(exists=True), required=True)
@click.option('--dm1label', type=str, required=True)
@click.option('--dm2label', type=str, required=True)
@click.option('--output', type=click.Path(exists=False), required=True)
def scatter(dm1, dm2, dm1label, dm2label, output):
    effdm1 = qiime2.Artifact.load(dm1).view(pd.DataFrame)
    effdm2 = qiime2.Artifact.load(dm2).view(pd.DataFrame)
    effdm1.set_index('column', inplace=True)
    effdm2.set_index('column', inplace=True)

    effdm1label = 'Effect size %s' % dm1label.replace('-', ' ')
    effdm2label = 'Effect size %s' % dm2label.replace('-', ' ')
    effdm1.rename(columns={'effect_size': effdm1label}, inplace=True)
    effdm2.rename(columns={'effect_size': effdm2label}, inplace=True)
    eff = effdm1.copy()
    eff[effdm2label] = effdm2[effdm2label]

    r, p = ss.pearsonr(eff[effdm1label], eff[effdm2label])

    r2 = r**2
    plt.figure(figsize=(10,10))
    sns.set(font_scale=2)
    sns.set_style(style='white')
    sns.scatterplot(data=eff, x=effdm1label, y=effdm2label, s=250)
    ax = plt.gca()

    min_ = 0
    max_ = max(eff[effdm1label].max(), eff[effdm2label].max()) + .02
    ax.set_xlim(min_, max_)
    ax.set_ylim(min_, max_)
    ax.plot([0, max_], [0, max_], color='r')

    ax.text(0.01, max_ - 0.02, 'Pearson r^2=%0.2f; p=%0.2e' % (r2, p), fontsize=24)
    ax.set_title("Beta-diversity effect sizes (Cohen's D and F)")

    if 'thdmi_cohort' in eff.index:
        focus_a = eff.loc['thdmi_cohort']
        x = focus_a[effdm1label]
        y = focus_a[effdm2label]
        ax.scatter([x], [y], marker='*', s=1000, c='g')

        focus_b = eff.loc['alcohol_types_white_wine']
        x = focus_b[effdm1label]
        y = focus_b[effdm2label]
        ax.scatter([x], [y], marker='*', s=1000, c='g')

        focus_c = eff.loc['alcohol_types_red_wine']
        x = focus_c[effdm1label]
        y = focus_c[effdm2label]
        ax.scatter([x], [y], marker='*', s=1000, c='g')
    plt.tight_layout()
    plt.savefig(output + '.pdf')


if __name__ == '__main__':
    scatter()
