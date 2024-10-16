rule miniprot_index:
    input:
        fasta = multiext('BPC/data/freeze_1/{sample}.fa.gz','','.fai','.gzi')
    output:
        index = 'analyses/pangene/{sample}.mpi'
    threads: 2
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '1h'
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
        runtime = '4h'
    shell:
        '''
miniprot -t {threads} --outs=0.97 -I -u {input.fasta} {input.peptides} |\
pigz -p {threads} -c > {output.paf}
        '''

rule pangene:
    input:
        paf = expand(rules.miniprot.output['gff'],sample=samples)
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
