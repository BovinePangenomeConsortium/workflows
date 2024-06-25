samples = glob_wildcards('assemblies/{sample}.fasta.gz').sample
config['reference'] = '/cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa'

rule ragtag_scaffold:
    input:
        fasta = multiext('assemblies/{sample}.fasta.gz','','.fai','.gzi'),
        reference = config['reference']
    output:
        fasta = multiext('pangenome/rescaffolded/{sample}.fasta.gz','','.fai','.gzi')
    params:
        mm2_opt = '-x asm20'
    threads: 8
    resources:
        mem_mb = 4000,
        scratch = '10G'
    shell:
        '''
        seqtk cutN -n 0 {input.fasta[0]} > $TMPDIR/asm.fa
        ragtag.py scaffold {input.reference} $TMPDIR/asm.fa -o $TMPDIR --mm2-params "{params.mm2_opt} -t {threads}" 

        sed 's/_RagTag//g' $TMPDIR/ragtag.scaffold.fasta |
        bgzip -@ {threads} -c > {output.fasta[0]}
        samtools faidx {output.fasta[0]}
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

