# Uses the layout defined in the BPC spreadsheet
def map_ID_to_filename(**wildcards):
    filename_map = dict(pl.read_csv('metdata.csv')
                       .select(['ID','Filename'])
                       .iter_rows())
    return filename_map[wildcards.sample]

rule panSN_renaming:
    input:
        lambda wildcards: expand(rules.ragtag_scaffold.output['fasta'][0],sample=map_ID_to_filename(wildcards.sample))
    output:
        'data/freeze_1/{sample}.fa.gz'
    shell:
        '''
        # we want to rename the chromosomes in panSN
        # we also want to rename the output file to match the ID schema
        # rename -> can we just modify ragtag based on the agp?
        '''

rule agc_create:
    input:
        assemblies = expand(rules.panSN_renaming.output['fasta'],samples=samples)
    output:
        agc = 'data/freeze_1/BPC_freeze_1.agc'
    threads: 8
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '4h'
    shell:
        '''
        agc create -d -t {threads} {input.assemblies} > {output.agc}
        '''

rule panSN_spec:
    input:
        expand(rules.ragtag_scaffold.output['fasta'][0],sample=samples)
    output:
        multiext('pangenome/input/{chromosome}.fa.gz','','.fai','.gzi')
    threads: 4
    resources:
        mem_mb = 1500,
        walltime = '30m'
    shell:
        '''
        for i in {input}
        do
          S=$(basename ${{i%.fasta.gz}})
          samtools faidx --length 0 $i {wildcards.chromosome} | fastix -p "$S#0#" -
        done | bgzip -@ {threads} -c > {output[0]}
        samtools faidx {output[0]}
        '''

rule mash_triangle:
    input:
        rules.panSN_spec.output
    output:
        'pangenome/tree/{chromosome}.matrix'
    threads: 4
    resources:
        mem_mb = 2500,
        walltime = '1h'
    shell:
        '''
        mash triangle -i -p {threads} {input[0]} > {output}
        '''
