config['reference'] = '/cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa'

import polars as pl
from pathlib import Path

def get_samples():
    metadata = pl.read_csv('summary.csv')#'£samplesheet_20240628.csv')
    return (metadata#.filter(pl.col('Species')=='Bos taurus')
                    .get_column('sample').to_list()
           )

samples = get_samples()

rule all:
    input:
        expand('analyses/scaffolding/{reference}/{sample}',reference=config['references'],sample=samples)

## We start by downloading the frozen set of assemblies
rule get_BPC_agc:
    output:
        metadata = expand('data/downloads/metadata.csv'),
        agc = expand('data/downloads/BPC.agc')
    localrule: True
    shell:
        '''
        wget
        '''

rule extract_BPC_agc:
    input:
        rules.get_BPC_agc.output
    output:
        expand('data/raw_assemblies/{sample}.fa',sample=samples)
    shell:
        '''
        agc extract
        '''

rule cut_assemblies_at_gaps:
    input:
        fasta = multiext('data/raw_assemblies/{sample}.fasta.gz','','.fai','.gzi')
    output:
        fasta = 'data/contigs/{sample}.fa'
    threads: 1
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '30m'
    shell:
        '''
        seqtk cutN -n 0 {input.fasta[0]} > {output}
        '''

rule ragtag_scaffold:
    input:
        fasta = rules.cut_assemblies_at_gaps.output['fasta'],
        reference = lambda wildcards: config['references'][wildcards.reference]
    output:
        multiext('analyses/scaffolding/{reference}/{sample}/ragtag.scaffold','.agp','.fasta','.err','.confidence.txt','.stats','.asm.paf','.asm.paf.log')
    params:
        _dir = lambda wildcards, output: Path(output[0]).parent,
        mm2_opt = '-x asm10'
    threads: 6
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '2h'
    shell:
        '''
        ragtag.py scaffold {input.reference} {input.fasta} -o {params._dir} -e {input.unplaced_contigs} --mm2-params "{params.mm2_opt} -t {threads}"
        '''

## compare agp files for multiple reference choices
