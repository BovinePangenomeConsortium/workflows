#subset by graph
def get_reference_sense_path():
    for entry in metadata.filter(pl.col('Reference annotation')=='Y').iter_rows(named=True):
        yield f"{entry['Animal ID']}"
        # see https://github.com/vgteam/vg/issues/4533
        #yield f"{entry['Animal ID']}#{entry['Haplotype']}"


#vg simplify -t 2 -P 85D7A68E -L 0.8 -m 20 small.pg > small_filtered_no_K.pg
## simplify graph upstream?

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
        mem_mb_per_cpu = 60000,
        runtime = '24h'
    shell:
        '''
echo {input.og} | tr ' ' '\\n' > $TMPDIR/og.fofn
odgi squeeze --optimize --threads {threads} --input-graphs $TMPDIR/og.fofn --out /dev/stdout |
tee {output.og} |
odgi view --threads {threads} --idx /dev/stdin --to-gfa > {output.gfa}
        '''

rule vg_combine:
    input:
        gfa = expand('analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.unchop.gfa',chromosome=ALL_CHROMOSOME,allow_missing=True)
    output:
        vg = 'analyses/pggb/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.vg'
    threads: 1
    resources:
        mem_mb_per_cpu = 60000,
        runtime = '4h'
    shell:
        '''
vg combine {input.gfa} > {output.vg}
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
        gbz = multiext('analyses/giraffe/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop','.dist','.giraffe.gbz','.longread.zipcodes','.longread.withzip.min','.shortread.zipcodes','.shortread.withzip.min')
    params:
        prefix = lambda wildcards, output: PurePath(output.gbz[0]).with_suffix('')
    threads: 1
    resources:
        mem_mb_per_cpu = 150000,
        runtime = '24h'
    shell:
        '''
vg autoindex --prefix {params.prefix} \
--threads {threads} \
--tmp-dir $TMPDIR \
--ref-fasta {input.fasta[0]} \
--gfa {input.gfa} \
--workflow sr-giraffe \
--workflow lr-giraffe
        '''

rule kmc_count:
    input:
        fastq = expand('/cluster/work/pausch/alex/genome_alignment/fastq/{sample}.hap1.R{N}.SR.fq.gz',N=(1,2),allow_missing=True)
        #fastq = 'data/sequencing_reads/{sample}.fq'
    output:
        kff = 'analyses/giraffe/kmers/{sample}.kff'
    params:
        prefix = lambda wildcards, output: PurePath(output['kff']).with_suffix('')
    threads: 4
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '1h'
    shell:
        '''
kmc -k29 -m32 -okff -t{threads} -hp -fq @<(echo {input.fastq} | tr " " "\\n") {params.prefix} $TMPDIR
        '''

rule vg_gbwt:
    input:
        #gfa = rules.odgi_squeeze.output['gfa'],
        vg = rules.vg_combine.output['vg']
    output:
        gbz = 'analyses/giraffe/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.direct.gbz',
        ri = 'analyses/giraffe/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.direct.ri'
    params:
        graph_in = lambda wildcards, input: "--gfa-input {input[0]}" if Path(input[0]).suffix == '.gfa' else f"--xg-name {input[0]}"
    threads: 2
    resources:
        mem_mb_per_cpu = 50000,
        runtime = '4h'
    shell:
        '''
vg gbwt \
{params.graph_in} \
--gbz-format \
--graph-name {output.gbz} \
--r-index {output.ri} \
--num-threads {threads}
        '''

rule vg_index:
    input:
        gbz = rules.vg_gbwt.output['gbz']
    output:
        dist = 'analyses/giraffe/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.direct.dist'
    threads: 2
    resources:
        mem_mb_per_cpu = 25000,
        runtime = '4h'
    shell:
        '''
vg index \
--temp-dir $TMPDIR \
--no-nested-distance \
--dist-name {output.dist} \
--threads {threads} \
{input.gbz}
        '''

rule vg_haplotype:
    input:
        gbz = rules.vg_gbwt.output['gbz'],
        ri = rules.vg_gbwt.output['ri'],
        dist = rules.vg_index.output['dist']
    output:
        hapl = 'analyses/giraffe/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.direct.hapl'
    params:
        reference_sense = ' '.join([f"--set-reference {path}" for path in get_reference_sense_path()])
    threads: 4
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '4h'
    shell:
        '''
vg haplotypes \
--diploid-sampling \
--include-reference \
{params.reference_sense} \
--threads {threads} \
--haplotype-output {output.hapl} \
{input.gbz}
        '''

rule vg_autoindex_personalised:
    input:
        hapl = '',
        gbz = '',
        kff = ''
    output:
        'something'
    shell:
        '''
vg haplotypes \
--threads {threads} \
--include-reference \
--diploid-sampling \
-i graph.hapl \
-k ${TMPDIR}/sample.kff \
-g ${TMPDIR}/sampled.gbz \
graph.gbz
        '''

#integrated approach
#--set-reference?
#reference fasta?
rule vg_giraffe_personalised:
    input:
        kff = rules.kmc_count.output['kff'],
        gbz = rules.vg_gbwt.output['gbz'],
        hapl = rules.vg_haplotype.output['hapl'],
        fastq = expand('/cluster/work/pausch/alex/genome_alignment/fastq/{sample}.hap1.R{N}.SR.fq.gz',N=(1,2),allow_missing=True)
    output:
        gam = 'analyses/giraffe/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.{sample}.personalised.gam'
    params:
        fastq = lambda wildcards, input: ' '.join([f'--fastq-in {fq}' for fq in input.fastq])
    threads: 8
    resources:
        mem_mb_per_cpu = 9000,
        runtime = '24h'
    shell:
        '''
vg giraffe \
--threads {threads} \
--gbz-name {input.gbz} \
--haplotype-name {input.hapl} \
--kff-name {input.kff} \
--parameter-preset default \
{params.fastq} > {output.gam}
        '''

rule fake:
    output:
        'fake'
    shell:
        '''
vg giraffe \
--haplotype-name {input.hapl} \
--kff-name {input.kff} \
--gbz-name {input.gbz} \
--index-basename {params.prefix} \
--sample {wildcards.sample} \
--include-reference {params.ref}
        '''
#vg giraffe -p -t 16 -Z graph.gbz --haplotype-name graph.hapl --kmer-name ${TMPDIR}/sample.kff -N sample -i -f sample.fq.gz > sample.gam

def get_sample_read_type(wildcards):
    match alignment_metadata.filter(pl.col('sample ID')==wildcards.sample).get_column('read type')[0]:
        case 'HiFi':
            return 'hifi'
        case 'Illumina-PE':
            return 'default'
        case 'Illumina-SE':
            return 'chaining-sr'
        case _:
            return 'default'

def get_sample_bam(wildcards):
    return alignment_metadata.filter(pl.col('sample ID')==wildcards.sample).get_column('bam')[0]

rule vg_giraffe_LR:
    input:
        gbz = rules.vg_autoindex.output['gbz'],
        bam = get_sample_bam
    output:
        gam = 'analyses/giraffe/{graph}/p{p}_s{segment_length}/whole_genome.k{k}.POA{POA}.unchop.{sample}.LR.gam'
    params:
        mode = get_sample_read_type
    threads: 8
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '4h'
    shell:
        '''
samtools fastq -@ {threads} {input.bam} |\
vg giraffe --parameter-preset {params.mode} --threads {threads} --gbz-name {input.gbz[1]} --fastq-in /dev/stdin > {output.gam}
        '''

rule vg_stats:
    input:
        gam = rules.vg_giraffe_LR.output['gam']
    output:
        stats = ''
    shell:
        '''
vg stats --alignments {input.gam} > {output.stats}
        '''
