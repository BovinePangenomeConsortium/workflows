import polars as pl
from pathlib import PurePath


ALL_CHROMOSOME = list(map(str,range(1,30))) + ['X','Y','MT']

wildcard_constraints:
    sample = r'[\w+\.\-_]+',
    graph = r'all'+ ''.join(['|'+graph for graph in config.get('pangenomes') if graph]),
    chromosome = r'|'.join(ALL_CHROMOSOME),
    L = r'\d+',
    reference = '|'.join(config.get('references')),
    peptides = '|'.join(config.get('peptides'))

metadata = pl.read_csv(config['metadata'],infer_schema_length=10000)
ANNOTATED_GENOMES = metadata.filter(pl.col('Reference annotation')=='Y').get_column('Filename').to_list()

ALL_CHROMOSOME = list(map(str,range(1,30))) + ['X','Y','MT']

alignment_metadata = pl.read_csv(config['alignment_metadata'],infer_schema_length=10000) if 'alignment_metadata' in config else pl.DataFrame()

def determine_pangenome_samples(wildcards):
    try:
        match wildcards.graph:
             case 'breed_representative':
                 subset = metadata.filter((pl.col('Species')=='Bos taurus')&(pl.col('Breed representative')=='Y'))
             case 'subspecies_representative':
                 subset = metadata.filter((pl.col('Species')=='Bos taurus')&(pl.col('Subspecies representative')=='Y'))
             case 'all':
                 subset = metadata
             case _:
                 subset = metadata.filter(pl.col('Filename').is_in(config['pangenomes'][wildcards.graph]))
    except AttributeError:
        subset = metadata 
    return subset.get_column('Filename').to_list()

include: 'snakepit/pangenome.smk'
include: 'snakepit/QC_assemblies.smk'
#include: 'snakepit/minigraph.smk'
#include: 'snakepit/pggb.smk'
#include: 'snakepit/pangene.smk'
#include: 'snakepit/pangenome_alignment.smk'
#include: 'snakepit/centomeres.smk'

rule all:
    input:
        'analyses/QC_summary.all.csv',
        #expand('analyses/pggb/medium_test/p95_s5000/{chromosome}.k31.POAasm20.unchop.graph.png',chromosome=ALL_CHROMOSOME)
