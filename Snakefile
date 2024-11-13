import polars as pl

include: 'snakepit/pangenome.smk'
include: 'snakepit/minigraph.smk'
include: 'snakepit/pggb.smk'

wildcard_constraints:
    graph = r'all|subspecies_representative|breed_representative'

def determine_pangenome_samples(graph):
    metadata = pl.read_csv(config['metadata'])

    match graph:
        case 'breed_representative':
            return (metadata.filter(pl.col('Species')=='Bos taurus')
                            .filter(pl.col('Breed representative')=='Y').get_column('ID').to_list()
                   )
        case 'subspecies_representative':
            return (metadata.filter(pl.col('Species')=='Bos taurus')
                            .filter(pl.col('Breed representative')=='Y').get_column('ID').to_list()
                   )
        case 'all' | _:
            return metadata.get_column('ID').to_list()

rule all:
    input:
        'analyses/minigraph/breed_representative/25.gfa',
        'analyses/pggb/breed_representative/p95_s5000/25.k23.POAasm5.unchop.gfa'
