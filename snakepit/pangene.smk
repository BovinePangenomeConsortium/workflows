rule miniprot_index:
    input:
        fasta = rules.panSN_renaming.output['fasta']
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
        paf = lambda wildcards: expand(rules.miniprot_align.output['paf'],sample=determine_pangenome_samples(wildcards.graph))
    output:
        gfa = 'analyses/pangene/{graph}.{reference}.gfa'
    threads: 1
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '1h'
    shell:
        '''
pangene {input.paf} > {output.gfa} 
        '''

rule pangene_matrix:
    input:
        gfa = rules.pangene.output['gfa']
    output:
        tsv = 'analyses/pangene/{graph}.{reference}.tsv'
    localrule: True
    shell:
        '''
pangene.js gfa2matrix -c {input.gfa} > {output.tsv}
        '''

rule pangene_call:
    input:
        gfa = rules.pangene.output['gfa']
    output:
        'analyses/pangene/{graph}.{reference}.call'
    localrule: True
    shell:
        '''
pangene.js call -p -s {input.gfa} > {output}
        '''
