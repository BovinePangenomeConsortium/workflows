import polars as pl

rule all:
    input:
        'analyses/minigraph/mash.txt'

rule mash:
    input:
        fasta = expand('data/freeze_1/{sample}.fa.gz',sample=pl.read_csv(config['metadata']).get_column('ID').to_list())
    output:
        'analyses/minigraph/mash.txt'
    threads: 4
    resources:
        mem_per_cpu_mb = 2500
    shell:
        '''
mash dist -p {threads} {input} > {output}
        '''
