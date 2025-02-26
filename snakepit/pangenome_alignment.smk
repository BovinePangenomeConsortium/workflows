rule odgi_squeeze:
    input:
        og = expand('analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.unchop.og',chromosome=ALL_CHROMOSOME,allow_missing=True)
    output:
        og = 'analyses/pggb/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.og',
        gfa = 'analyses/pggb/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.gfa'
    threads: 2
    resources:
        mem_mb_per_cpu = 20000,
        walltime = '4h'
    shell:
        '''
echo {input.og} | tr ' ' '\\n' > $TMPDIR/og.fofn
odgi squeeze --optimize -t {threads} -f $TMPDIR/og.fofn -o /dev/stdout |
tee {output.og} |
odgi view -i /dev/stdin -g > {output.gfa} #TODO: sed the reference header here
        '''

rule concat_reference_sequence:
    input:
        fasta = expand(rules.panSN_gather.output['fasta'][0],chromosome=ALL_CHROMOSOME,allow_missing=True)
    output:
        fasta = multiext('data/currated_assemblies/{graph}/whole_genome.fa.gz','','.fai','.gzi')
    localrule: True
    shell:
        '''
cat {input} > {output[0]}
samtools faidx {output[0]}
        '''

rule vg_autoindex:
    input:
        gfa = rules.odgi_squeeze.output['gfa'],
        fasta = rules.concat_reference_sequence.output['fasta']
    output:
        gbz = multiext('analyses/pggb/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop','.gbz','.dist')
    params:
        prefix = lambda wildcards, output: PurePath(output.gbz[0]).with_suffix('')
    threads: 1
    resources:
        mem_mb_per_cpu = 200000,
        walltime = '24h'
    shell:
        '''
vg autoindex -p {params.prefix} -t {threads} -T $TMPDIR -r {input.fasta[0]} -g {input.gfa} -w lr-giraffe
        '''

rule kmc_count:
    input:
        fastq = 'data/sequencing_reads/{sample}.fq'
    output:
        kff = multiext('analyses/pangenome/giraffe/{sample}','.kmc_pre','.kmc_suf')
    params:
        prefix = lambda wildcards, output: PurePath(output['kff']).with_suffix('')
    shell:
        '''
kmc -k29 -m32 -okff -t {threads} -hp {input.fastq} {params.prefix} $TMPDIR
        '''

rule vg_giraffe_haplotype:
    input:
        kff = rules.kmc_count.output['kff'],
        hapl = '',
        gbz = ''
    output:
        ''
    shell:
        '''
vg giraffe --haplotype-name {input.hapl} --kff-name {input.kff} \
--gbz-name {input.gbz} --index-basename {params.prefix} \
--sample {wildcards.sample} --include-reference {params.ref}
        '''
#vg giraffe -p -t 16 -Z graph.gbz --haplotype-name graph.hapl --kmer-name ${TMPDIR}/sample.kff -N sample -i -f sample.fq.gz > sample.gam

rule vg_giraffe_LR:
    input:
        ''
    output:
        gam = ''
    shell:
        '''
vg giraffe -b hifi -Z hprc-v1.1-mc-chm13.d9.gbz -f longread/hifi.fq -p >hifi.mapped.gam
        '''

rule vg_haplotype:
    input:
        ''
    output:
        ''
    shell:
        '''
        vg index -j graph.dist --no-nested-distance graph.gbz
vg gbwt -p --num-threads 16 -r graph.ri -Z graph.gbz
vg haplotypes --diploid-sampling --include-reference -v 2 -t 16 -H graph.hapl graph.gbz
        '''

rule vg_stats:
    input:
        gam = rules.vg_giraffe_LR.output['gam']
    output:
        stats = ''
    shell:
        '''
vg stats -a {input.gam} > {output.stats}
        '''
