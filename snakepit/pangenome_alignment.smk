def get_reference_sense_path():
    for entry in metadata.filter(pl.col('Reference annotation')=='Y').iter_rows(named=True):
        yield f"{entry['Animal ID']}"
        # see https://github.com/vgteam/vg/issues/4533
        #yield f"{entry['Animal ID']}#{entry['Haplotype']}"

rule odgi_squeeze:
    input:
        og = expand('analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.unchop.og',chromosome=ALL_CHROMOSOME,allow_missing=True)
    output:
        og = 'analyses/pggb/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.og',
        gfa = 'analyses/pggb/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.gfa'
    params:
        reference_IDs = f"RS:Z:{' '.join(get_reference_sense_path())}"
    threads: 2
    resources:
        mem_mb_per_cpu = 20000,
        runtime = '4h'
    shell:
        '''
echo {input.og} | tr ' ' '\\n' > $TMPDIR/og.fofn
odgi squeeze --optimize -t {threads} -f $TMPDIR/og.fofn -o /dev/stdout |
tee {output.og} |
odgi view -i /dev/stdin -g |\
sed '1s/$/\\t{params.reference_IDs}/' > {output.gfa}
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
        gbz = multiext('analyses/pggb/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop','.dist','.giraffe.gbz','.longread.zipcodes','.longread.withzip.min','.shortread.zipcodes','.shortread.withzip.min')
    params:
        prefix = lambda wildcards, output: PurePath(output.gbz[0]).with_suffix('')
    threads: 1
    resources:
        mem_mb_per_cpu = 150000,
        runtime = '24h'
    shell:
        '''
vg autoindex -p {params.prefix} -t {threads} -T $TMPDIR -r {input.fasta[0]} -g {input.gfa} -w sr-giraffe -w lr-giraffe
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
        gbz = rules.vg_autoindex.output['gbz'],
        bam = ''
    output:
        gam = 'analyses/giraffe/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.{sample}.gam'
    threads: 8
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '4h'
    shell:
        '''
samtools fastq -@ {threads} {input.bam} |\
vg giraffe --parameter-preset hifi --threads {threads} --gbz-name {input.gbz[1]} --fastq-in /dev/stdin > {output.gam}
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
