rule miniprot_index:
    input:
        config['reference']
    output:
        index = 'analyses/pangene/{reference}.mpi'
    threads: 4
    shell:
        '''
miniprot -t {threads} -d {output.index} {input.fasta}
        '''

rule miniprot_align:
    input:
        fasta = multiext('BPC/data/freeze_1/{sample}.fa.gz','','.fai','.gzi'),
        cDNAs = lambda wildcards: config['cDNAs'][wildcards.reference]
    output:
        paf = 'analyses/pangene/{sample}.{reference}.paf.gz'
    threads: 8
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '4h'
    shell:
        '''
miniprot -t {threads} {input.fasta[0]} {input.cDNAs} |\
pigz -p {threads} -c > {output.paf}
        '''

rule pangene:
    input:
        paf = expand(rules.miniprot.output['gff'],sample=samples)
    output:
        gfa = 'analyses/pangene/{reference}.gfa'
    shell:
        '''
pangene {input.paf} > {output.gfa} 
        '''
