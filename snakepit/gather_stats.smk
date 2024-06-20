samples = glob_wildcards('assemblies/{sample}.fasta.gz').sample
config['reference]'] = '/cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa'

rule all:
    input:
        'summary.csv'

rule calculate_N50:
    input:
        fasta = multiext('assemblies/{sample}.fasta.gz','','.fai','.gzi')
    output:
        'contiguity/{sample}.N50'
    resources:
        walltime = '10m'
    shell:
        '''
        seqtk cutN -n 0 {input.fasta[0]} | calN50.js -L 2770686122 - > {output}
        '''

rule calculate_gene_completeness:
    input:
        fasta = multiext('assemblies/{sample}.fasta.gz','','.fai','.gzi')
    output:
        _dir = directory('completeness/compleasm_{sample}'),
        full_table = 'completeness/compleasm_{sample}/cetartiodactyla_odb10/full_table.tsv'
    threads: 4
    resources:
        mem_mb = 12500,
        walltime = '2h'
    shell:
        '''
        compleasm run -a {input.fasta[0]} -o {output._dir} -l cetartiodactyla -t {threads}
        '''

rule minimap2_reference_aligned:
    input:
        fasta = multiext('assemblies/{sample}.fasta.gz','','.fai','.gzi')
    output:
        'reference_alignment/{sample}.ARS_UCD2.0.paf'
    threads: 4
    resources:
        mem_mb = 10000,
        walltime = '2h'
    shell:
        '''
        minimap2 -t {threads} --cs -cxasm10 {config[reference]} {input.fasta[0]} > {output}
        '''

rule calculate_variant_level:
    input:
        paf = rules.minimap2_reference_aligned.output
    output:
        vcf = multiext('reference_alignment/{sample}.vcf.gz','','.tbi')
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '1h'
    shell:
        '''
        sort -k6,6 -k8,8n {input.paf} |\
        paftools call -f {config[reference]} -s {wildcards.sample} - |\
        bcftools view --write-index -o {output.vcf[0]}
        '''

rule summarise_sample_metrics:
    input:
        N50 = rules.calculate_N50.output,
        comepleteness = rules.calculate_gene_completeness.output['full_table'],
        variants = rules.calculate_variant_level.output['vcf'],
    output:
        'summary/{sample}.csv'
    shell:
        '''
        echo -n "{wildcards.sample}," > {output}
        awk '$1~/(SZ|NN)/ {{printf $2",";next}} {{if ($1=="NL"&&$2==50) {{printf $3","} }}}' {input.N50}' >> {output}

        '''

rule summarise_all_metrics:
    input:
        expand(rules.summarise_sample_metrics.output,sample=samples)
    output:
        'summary.csv'
    localrule: True
    shell:
        '''
        cat {input} > {output}
        '''
