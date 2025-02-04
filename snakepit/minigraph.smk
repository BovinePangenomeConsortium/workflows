rule mash_sketch:
    input:
        fasta = multiext('data/currated_assemblies/{sample}.fa.gz','','.fai','.gzi')
    output:
        sketch = 'analyses/minigraph/mash/{sample}.msh'
    params:
        prefix = lambda wildcards, output: PurePath(output['sketch']).with_suffix('')
    threads: 2
    resources:
        mem_per_cpu_mb = 2500,
        runtime = '1h'
    shell:
        '''
mash sketch -p {threads} -o {params.prefix} {input.fasta[0]}
        '''

rule mash_dist:
    input:
        lambda wildcards: expand(rules.mash_sketch.output['sketch'],sample=determine_pangenome_samples(wildcards.graph))
    output:
        'analyses/minigraph/{graph}.mash.txt'
    threads: 4
    resources:
        mem_per_cpu_mb = 2500,
        runtime = '1h'
    shell:
        '''
mash dist -p {threads} {input} > {output}
        '''

rule minigraph_construct:
    input:
        assemblies = lambda wildcards: expand('data/currated_assemblies/chromosomes/{sample}.{chromosome}.fa.gz',sample=determine_pangenome_samples(wildcards.graph),allow_missing=True)
    output:
        gfa = 'analyses/minigraph/{graph}/L{L}/{chromosome}.basic.gfa'
    threads: 1
    resources:
        mem_mb_per_cpu_per_cpu = 20000,
        runtime = '24h'
    shell:
        '''
minigraph -t {threads} -cxggs -j 0.05 -L {wildcards.L} {input.assemblies} > {output.gfa}
        '''

rule minigraph_call:
    input:
        gfa = rules.minigraph_construct.output['gfa'],
        assembly = 'data/currated_assemblies/chromosomes/{sample}.{chromosome}.fa.gz'
    output:
        bed = 'analyses/minigraph/{graph}/L{L}/{sample}.{chromosome}.bed'
    threads: 1
    resources:
        mem_mb_per_cpu_per_cpu = 10000,
        runtime = '1h'
    shell:
        '''
minigraph -t {threads} -cxasm --call -j 0.05 -L {wildcards.L} {input.gfa} {input.assembly} > {output.bed}
        '''

rule minigraph_path:
    input:
        paths = lambda wildcards: expand(rules.minigraph_call.output['bed'],sample=determine_pangenome_samples(wildcards.graph),allow_missing=True),
        gfa = rules.minigraph_construct.output['gfa']
    output:
        gfa = 'analyses/minigraph/{graph}/L{L}/{chromosome}.gfa'
    threads: 1
    resources:
        mem_mb_per_cpu_per_cpu = 5000,
        runtime = '1h'
    shell:
        '''
#needs special branch of mgutils.js
{{ cat {input.gfa} ; paste {input.paths} | mgutils.js path - ; }} > {output}
        '''
