rule KMC_count:
    input:
        fasta = rules.panSN_renaming.output['fasta']#rules.samtools_fastq.output,
        #coverage = rules.estimate_coverage.output
    output:
        kmers = multiext('analyses/satellites/{sample}.kmc','.kmc_pre','.kmc_suf','.txt')
    params:
        prefix = lambda wildcards, output: Path(output[0]).with_suffix(''),
        #threshold = lambda wildcards, input: int(float(open(input.coverage[0]).read())*10)
    threads: 2
    resources:
        mem_mb_per_cpu = 10000,
        walltime = '1h'
    shell:
        '''
kmc -k151 -t{threads} -ci20 -cs100000 -fm {input.fasta[0]} {params.prefix} $TMPDIR

kmc_dump {params.prefix} {output.kmers[2]}
        '''

rule SRF:
    input:
        kmers = rules.KMC_count.output['kmers'][2]
    output:
        satellites = 'analyses/satellites/{sample}.srf.fa'
    resources:
        mem_mb_per_cpu = 2500,
        walltime = '30m'
    shell:
        '''
        srf -p {wildcards.sample} -l 40 -c 100 {input.kmers} > {output.satellites}
        '''

rule srf_gather:
    input:
        paf = expand(rules.SRF.output['satellites'],sample=determine_pangenome_samples,allow_missing=True)
    output:
        gfa = 'analyses/satellites/{graph}.txt'
    localrule: True
    shell:
        '''
        touch {output}
        '''
