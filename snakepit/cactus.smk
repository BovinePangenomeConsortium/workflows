# very stripped back version of https://github.com/harvardinformatics/cactus-snakemake/blob/main/cactus.smk

#TODO: optimize resource usage and add more config options
#TODO: check on split and gather steps

rule all:
    input:
        'cactus/all/final.gfa.gz'

rule cactus_minigraph:
    input:
        seqfile = config['seq_file']
    output:
        gfa = 'cactus/{graph}/sv.gfa'
    threads: 8
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '4h'
    container: config["cactus_image"]
    shell: '''
cactus-minigraph $TMPDIR/minigraph {input.seqfile} {output.gfa} --reference {config[ref_genome]}
'''

rule cactus_graphmap:
    input:
        seqfile = config['seq_file'],
        gfa = rules.cactus_minigraph.output.gfa
    output:
        multiext('cactus/{graph}/graphmap', paf='.paf', fasta='.sv.gfa.fa', gaf='.gaf.gz', paf_filter_log='.paf.filter.log', paf_unfiltered='.paf.unfiltered.gz')
    threads: 8
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '4h'
    container: config["cactus_image"]
    shell: '''
cactus-graphmap $TMPDIR/graphmap {input.seqfile} {input.gfa} {output.paf} --outputFasta {output.fasta} --reference {config[ref_genome]}
'''

checkpoint cactus_split:
    input:
        seqfile = config['seq_file'],
        gfa = rules.cactus_minigraph.output.gfa,
        paf = rules.cactus_graphmap.output.paf
    output:
        split_dir = directory( 'cactus/{graph}/split'),
    threads: 8
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '4h'
    container: config["cactus_image"]
    shell: '''
cactus-graphmap-split $TMPDIR/split {input.seqfile} {input.gfa} {input.paf} --outDir {output.split_dir} --reference {config[ref_genome]}
'''

rule cactus_align:
    input:
        seqfile = 'cactus/{graph}/split/{chromosome}.seqfile',
        paf = 'cactus/{graph}/split/{chromosome}.paf'
    output:
        multiext('cactus/{graph}/align/{chromosome}', hal='.hal', vg='.vg')
    threads: 8
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '4h'
    container: config["cactus_image"]
    shell: '''
cactus-align $TMPDIR/align {input.seqfile} {input.paf} {output.hal} --pangenome --reference {config[ref_genome]} --outVG
'''

def gatherChromosomes(wildcards):

    with open(checkpoints.cactus_split.get(**wildcards).output.contig_sizes) as f:
        return [line.strip().split("\t")[0] for line in f if line.endswith('avg')]

rule cactus_join:
    input:
        hal = expand(rules.cactus_align.output.hal, chromosome=gatherChromosomes),
        vg = expand(rules.cactus_align.output.vg, chromosome=gatherChromosomes)
    output:
        multiext('cactus/{graph}/final', hal='.full.hal', gfa='.gfa.gz', vcf='.vcf.gz', vcf_index='.vcf.gz.tbi', gbz='.gbz', raw_vcf='.raw.vcf.gz', raw_vcf_index='.raw.vcf.gz.tbi', stats='.stats.tgz')
    threads: 8
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '4h'
    container: config["cactus_image"]
    shell: '''
cactus-graphmap-join $TMPDIR/join --vg {input.vg} --hal {input.hal} --outDir split --outName test --reference {config[ref_genome]} --vcf --giraffe full clip
'''