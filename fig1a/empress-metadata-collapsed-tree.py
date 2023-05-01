import click
import pandas as pd


@click.command()
@click.option('--collapse-map', type=click.Path(exists=True), required=True)
@click.option('--tip-association', type=click.Path(exists=True), required=True)
@click.option('--output', type=click.Path(exists=False), required=True)
def process(collapse_map, tip_association, output):
    collapse_map = pd.read_csv(collapse_map, sep='\t').set_index('feature-id')
    tip_association = pd.read_csv(tip_association, sep='\t').set_index('feature-id')

    metadata = []
    for node, nodegrp in collapse_map.groupby('collapsed-tip'):
        try:
            observed_tips = tip_association.loc[list(nodegrp.index)]
        except:
            print(node, len(nodegrp))
            print(nodegrp.index)
            raise
        is_agp = observed_tips['is_agp'].any()
        is_emp = observed_tips['is_emp'].any()
        both = is_agp & is_emp
        phy = observed_tips['phylum_top20'].value_counts().index[0]

        if both:
            label = 'Both'
        elif is_agp:
            label = 'AGP'
        elif is_emp:
            label = 'EMP'
        else:
            label = 'Neither'

        metadata.append([node, is_agp, is_emp, label, phy])

    metadata = pd.DataFrame(metadata, columns=['feature-id', 'is_agp', 'is_emp', 'label', 'phylum_top20'])
    metadata.to_csv(output, sep='\t', index=False, header=True)

if __name__ == '__main__':
    process()
