
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

rule pggb_construct:
    input:
        fasta = rules.panSN_spec.output
    output:
        graphs = multiext('pangenome/graphs/{chromosome}/pggb','.gfa','.og')
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '4h',
        scratch = '50G'
    params:
        _dir = lambda wildcards, output: Path(output.graphs[0]).parent,
        divergence = 95,#lambda wildcards, input: read_mash_triangle(input.mash[0],True),
        min_match = 23, #pggb default
        segment_length = '25k'
    shell:
        '''
        pggb -i {input.fasta[0]} -o {params._dir} -t {threads} \
        -s {params.segment_length} -p {params.divergence} -k {params.min_match} \
        --skip-viz --temp-dir $TMPDIR

        mv {params._dir}/{wildcards.chromosome}.*.smooth.final.gfa {output.graphs[0]}
        mv {params._dir}/{wildcards.chromosome}.*.smooth.final.og {output.graphs[1]}
        '''
