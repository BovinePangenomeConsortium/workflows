import polars as pl

include: 'snakepit/pangenome.smk'
include: 'snakepit/QC_assemblies.smk'
include: 'snakepit/minigraph.smk'
include: 'snakepit/pggb.smk'
include: 'snakepit/pangene.smk'

wildcard_constraints:
    graph = r'all|subspecies_representative|breed_representative',
    #reference = '|'.join(config.get('peptides'))

metadata = pl.read_csv(config['metadata'])
ANNOTATED_GENOMES = metadata.filter(pl.col('Reference annotation')=='Y').get_column('Filename').to_list()

def determine_pangenome_samples(graph=None):

    match graph:
        case 'breed_representative':
            return (metadata.filter(pl.col('Species')=='Bos taurus')
                            .filter(pl.col('Breed representative')=='Y').get_column('Filename').to_list()
                   )
        case 'subspecies_representative':
            return (metadata.filter(pl.col('Species')=='Bos taurus')
                            .filter(pl.col('Breed representative')=='Y').get_column('Filename').to_list()
                   )
        case 'all' | _:
            return metadata.get_column('Filename').to_list()

rule all:
    input:
        'analyses/freeze_1/QC_summary.csv',
        'analyses/minigraph/breed_representative/25.gfa',
        'analyses/pggb/breed_representative/p95_s5000/25.k23.POAasm5.unchop.gfa'
