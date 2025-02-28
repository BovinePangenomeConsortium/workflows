import polars as pl
from pathlib import PurePath

wildcard_constraints:
    sample = r'[\w+\.\-_]+',
    graph = r'small_test|all|subspecies_representative|breed_representative',
    chromosome = r'\d+|X|Y|MT',
    L = r'\d+',
    reference = '|'.join(config.get('references'))

metadata = pl.read_csv(config['metadata'],infer_schema_length=10000)
ANNOTATED_GENOMES = metadata.filter(pl.col('Reference annotation')=='Y').get_column('Filename').to_list()

ALL_CHROMOSOME = list(map(str,range(1,30))) + ['X','Y','MT']

def determine_pangenome_samples(graph=None):
    subset = None
    match graph:
        case 'small_test':
            subset = metadata.filter(pl.col('Filename').is_in(['ARS-UCD1.3','CxR_raft_trioUL.dip.charolais-sire.p_ctg','Charolais.haplotype1.chrNames.20231003','Wagyu_haplotype1_v1.polished','ASM4388211v1']))
        case 'breed_representative':
            subset = (metadata.filter(pl.col('Species')=='Bos taurus')
                            .filter(pl.col('Breed representative')=='Y'))
        case 'subspecies_representative':
            subset = (metadata.filter(pl.col('Species')=='Bos taurus')
                            .filter(pl.col('Breed representative')=='Y'))
        case 'all' | _:
            subset = metadata
    return subset.get_column('Filename').to_list()

include: 'snakepit/pangenome.smk'
include: 'snakepit/QC_assemblies.smk'
include: 'snakepit/minigraph.smk'
include: 'snakepit/pggb.smk'
include: 'snakepit/pangene.smk'
include: 'snakepit/pangenome_alignment.smk'

rule all:
    input:
        'analyses/freeze_1/QC_summary.csv',
        'analyses/minigraph/breed_representative/25.gfa',
        'analyses/pggb/breed_representative/p95_s5000/25.k23.POAasm5.unchop.gfa'
