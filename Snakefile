import polars as pl
from pathlib import PurePath


ALL_CHROMOSOME = list(map(str,range(1,30))) + ['X','Y','MT']

wildcard_constraints:
    sample = r'[\w+\.\-_]+',
    graph = r'|'.join(config.get('pangenomes','all')),
    chromosome = r'|'.join(ALL_CHROMOSOME),
    L = r'\d+',
    reference = '|'.join(config.get('references'))

metadata = pl.read_csv(config['metadata'],infer_schema_length=10000)
ANNOTATED_GENOMES = metadata.filter(pl.col('Reference annotation')=='Y').get_column('Filename').to_list()

ALL_CHROMOSOME = list(map(str,range(1,30))) + ['X','Y','MT']

def determine_pangenome_samples(wildcards):
    
    subset = metadata.filter(pl.col('Filename').is_in(config['pangenomes'][wildcards.graph])) if ('pangenomes' in config and wildcards.graph != 'all') else metadata
    return subset.get_column('Filename').to_list()
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
        'analyses/QC_summary.all.csv',
        expand('analyses/pggb/medium_test/p95_s5000/{chromosome}.k31.POAasm20.unchop.graph.png',chromosome=ALL_CHROMOSOME)
