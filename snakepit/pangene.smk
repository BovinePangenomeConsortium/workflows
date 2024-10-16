import polars as pl

wildcard_constraints:
    reference = '|'.join(config.get('peptides'))

rule all:
    input:
        'analyses/pangene/ARS_UCD1.3.gfa'

rule miniprot_index:
    input:
        fasta = multiext('data/freeze_1/{sample}.fa.gz','','.fai','.gzi')
    output:
        index = 'analyses/pangene/{sample}.mpi'
    threads: 4
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '30m'
    shell:
        '''
miniprot -t {threads} -d {output.index} {input.fasta}
        '''

rule download_annotated_peptides:
    output:
        fasta = 'analyses/pangene/peptides/{reference}.pep.fa.gz'
    params:
        url = lambda wildcards: config['peptides'][wildcards.reference]
    localrule: True
    shell:
        '''
wget -O {output.fasta} {params.url} 
        '''

rule miniprot_align:
    input:
        fasta = rules.miniprot_index.output['index'],
        peptides = rules.download_annotated_peptides.output['fasta'] 
    output:
        paf = 'analyses/pangene/{sample}.{reference}.paf.gz'
    threads: 8
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '2h'
    shell:
        '''
miniprot -t {threads} --outs=0.97 -I -u {input.fasta} {input.peptides} |\
pigz -p {threads} -c > {output.paf}
        '''

rule pangene:
    input:
        paf = expand(rules.miniprot_align.output['paf'],sample=pl.read_csv(config['metadata']).get_column('ID').to_list(),allow_missing=True)
    output:
        gfa = 'analyses/pangene/{reference}.gfa'
    threads: 1
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '4h'
    shell:
        '''
pangene {input.paf} > {output.gfa} 
        '''
