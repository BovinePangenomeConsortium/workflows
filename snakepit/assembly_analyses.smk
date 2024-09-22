config['reference'] = '/cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa'

import polars as pl

def get_samples():
    metadata = pl.read_csv('samplesheet_20240628.csv')
    return metadata.filter(pl.col('Species')=='Bos taurus').get_column('ID').to_list()

samples = get_samples()

## We start by downloading the frozen set of assemblies
rule get_BPC_agc:
    output:
        metadata = expand('data/downloads/metadata.csv')
        agc = expand('data/downloads/BPC.agc)
    localrule: True
    shell:
        '''
        wget
        '''

rule extract_BPC_agc:
    input:
        rules.get_BPC_agc
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
        mem_mb = 2500,
        walltime = '30m'
    shell:
        '''
        seqtk cutN -n 0 {input.fasta[0]} > {output}
        '''

rule ragtag_scaffold:
    input:
        fasta = rules.cut_assemblies_at_gaps.output['fasta']
        reference = lambda wildcards: config['references'][wildcards.reference]
    output:
        fasta = multiext('analyses/scaffolding/{sample}.{reference}.fasta.gz','','.fai','.gzi'),
        _dir = directory('analyses/scaffolding/{reference}/{sample}')
    params:
        mm2_opt = '-x asm10'
    threads: 6
    resources:
        mem_mb = 8000,
        walltime = '2h'
    shell:
        '''
        ragtag.py scaffold {input.reference} {input.fasta} -o {output._dir} --mm2-params "{params.mm2_opt} -t {threads}"

        sed 's/_RagTag//g' $TMPDIR/ragtag.scaffold.fasta |
        bgzip -@ {threads} -c > {output.fasta[0]}
        samtools faidx {output.fasta[0]}
        '''

## compare agp files for multiple reference choices
