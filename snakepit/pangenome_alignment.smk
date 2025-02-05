rule kmc_count:
    input:
        fastq = ''
    output:
        kff = multiext('analyses/pangenome/giraffe/{sample}','.kmc_pre','.kmc_suf')
    params:
        prefix = lambda wildcards, output: PurePath(output['kff']).with_suffix('')
    shell:
        '''
kmc -k29 -m32 -okff -t {threads} -hp {input.fastq} {params.prefix} $TMPDIR
        '''

rule vg_autoindex:
    input:
        gfa = ''
    output:
        gbz = ''
    shell:
        '''
    giraffe, sr-giraffe, lr-giraffe
        ''''

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
        ''
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
        gam = rules.vg_giraffe.output['gam']
    output:
        stats = ''
    shell:
        '''
vg stats -a {input.gam} > {output.stats}
        '''