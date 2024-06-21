samples = glob_wildcards('assemblies/{sample}.fasta.gz').sample
config['reference'] = '/cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa'

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
        fasta = multiext('assemblies/{sample}.fasta.gz','','.fai','.gzi'),
        reference = config['reference']
    output:
        'reference_alignment/{sample}.ARS_UCD2.0.paf'
    threads: 4
    resources:
        mem_mb = 10000,
        walltime = '2h'
    shell:
        '''
        minimap2 -t {threads} --cs -cxasm10 {input.reference} {input.fasta[0]} > {output}
        '''

rule calculate_variant_level:
    input:
        paf = rules.minimap2_reference_aligned.output,
        reference = config['reference']
    output:
        vcf = multiext('reference_alignment/{sample}.vcf.gz','','.csi')
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '1h'
    shell:
        '''
        sort -k6,6 -k8,8n {input.paf} |\
        paftools.js call -f {input.reference} -s {wildcards.sample} - |\
        bcftools view --write-index -o {output.vcf[0]}
        '''

rule calculate_reference_coverage:
    input:
        paf = rules.minimap2_reference_aligned.output,
        fai = config['reference'] + ".fai"
    output:
        bed = 'reference_alignment/{sample}.covered.bed'
    threads: 1
    resources:
        mem_mb = 2500,
        walltime = '30m'
    shell:
        '''
        cut -f 6,8,9 {input.paf} |\
        bedtools sort -faidx {input.fai} -i /dev/stdin |\
        bedtools merge -d 0 -i /dev/stdin |\
        bedtools genomecov -g {input.fai} -i /dev/stdin |\
        awk -v OFS='\\t' '$2==0 {{S[$1]=$4; U[$1]=$5; next}} {{C[$1]=$5}} END {{for (k in C) {{print k,"0",S[k],"{wildcards.sample}",C[k],U[k]}} }}' |\
        grep -P "^(\d|X|Y|MT)" |\
        bedtools sort -faidx {input.fai} -i /dev/stdin > {output.bed}
        '''

rule summarise_sample_metrics:
    input:
        N50 = rules.calculate_N50.output,
        comepleteness = rules.calculate_gene_completeness.output['full_table'],
        variants = rules.calculate_variant_level.output['vcf'],
        bed = rules.calculate_reference_coverage.output['bed']
    output:
        'summary/{sample}.csv'
    resources:
        walltime = '10m'
    shell:
        '''
        echo -n "{wildcards.sample}," > {output}
        awk '$1~/(SZ|NN)/ {{printf $2",";next}} {{if ($1=="NL"&&$2==50) {{printf $3","}} }}' {input.N50} >> {output}
    
        #busco analysis
        bcftools stats {input.variants[0]} | awk '$1=="SN"&&$5~/(SNPs|indels)/ {{printf $6","}}' >> {output}

        awk -v OFS=',' '$1~/[[:digit:]]/ {{A+=$5;++n;next}} {{B[$1]=$5}} END {{print A/n,B["X"],B["Y"]}}' {input.bed} >> {output}
        '''

rule summarise_all_metrics:
    input:
        expand(rules.summarise_sample_metrics.output,sample=samples)
    output:
        'summary.csv'
    localrule: True
    shell:
        '''
        echo "sample,genome size, N contigs, NG50, SNPs, InDels"
        cat {input} > {output}
        '''
