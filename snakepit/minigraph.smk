import polars as pl
from pathlib import PurePath

rule all:
    input:
        'analyses/minigraph/mash.txt'

rule mash_sketch:
    input:
        fasta = multiext('data/freeze_1/{sample}.fa.gz','','.fai','.gzi')
    output:
        sketch = 'analyses/minigraph/mash/{sample}.msh'
    params:
        prefix = lambda wildcards, output: PurePath(output['sketch']).with_suffix('')
    threads: 2
    resources:
        mem_per_cpu_mb = 2500,
        walltime = '1h'
    shell:
        '''
mash sketch -p {threads} -o {params.prefix} {input.fasta[0]}
        '''

rule mash_dist:
    input:
        expand(rules.mash_sketch.output['sketch'],sample=pl.read_csv(config['metadata']).get_column('ID').to_list())
    output:
        'analyses/minigraph/mash.txt'
    threads: 4
    resources:
        mem_per_cpu_mb = 2500,
        walltime = '1h'
    shell:
        '''
mash dist -p {threads} {input} > {output}
        '''

rule minigraph_construct:
    input:
        assemblies = expand('data/{{chromosome}}/{sample}.fa', sample=pangenome_samples)
    output:
        gfa = 'analyses/minigraph/{graph}/{chromosome}.basic.gfa'
    threads: 1
    resources:
        mem_mb_per_cpu = 20000,
        runtime = '24h'
    params:
        sample_order = lambda wildcards, input:make_minigraph_order(input.mash_distances[0],input.assemblies),
    shell:
        '''
minigraph -t {threads} -cxggs -j 0.05 -L 30 {params.sample_order} > {output.gfa}
        '''

rule minigraph_call:
    input:
        gfa = rules.minigraph_construct.output['gfa'],
        sample = 'assemblies/{chromosome}/{sample}.fa'
    output:
        bed = 'analyses/minigraph/{graph}/{sample}.{chromosome}.bed'
    threads: 1
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '1h'
    shell:
        '''
minigraph -t {threads} -cxasm --call -j 0.05 -L 30 {input.gfa} {input.sample} > {output.bed}
        '''

localrules: minigraph_path
rule minigraph_path:
    input:
        paths = expand(rules.minigraph_call.output['bed'],sample=samples,allow_missing=True)
        paths = expand('graphs/minigraph/{{chromosome}}.{sample}.bed',sample=filter(lambda x: x != get_reference_ID(), pangenome_samples)),
        gfa = rules.minigraph_construct.output['gfa']
    output:
        gfa = 'analyses/minigraph/{graph}/{chromosome}.gfa'
    params:
        samples = '\\n'.join(filter(lambda x: x != get_reference_ID(), pangenome_samples))
    threads: 1
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '1h'
    shell:
        '''
#needs special branch of mgutils.js
{{ vg convert -r 0 -g {input.gfa} -f ; paste {input.paths} | mgutils.js path <(echo -e "{params.samples}") - | sed 's/s//g' ; }} > {output}
        '''
