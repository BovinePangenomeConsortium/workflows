#many assemblies are gzipped, but we want bgzipped for better random access
rule rebgzip_assemblies:
    input:
        'downloaded_assemblies/{sample}.fa.gz'
    output:
        multiext('data/raw_assemblies/{sample}.fasta.gz','','.fai','.gzi')
    threads: 4
    resources:
        mem_mb_per_cpu = 2000,
        runtime = '1h'
    shell:
        '''
pigz -dc -p 4 {input} |\
bgzip -@ 4 -c > {output[0]}

samtools index -@ 4 {output[0]}
        '''

rule calculate_N50:
    input:
        fasta = rules.panSN_renaming.output['fasta']
    output:
        'analyses/QC/contiguity/{sample}.N50'
    resources:
        runtime = '10m'
    shell:
        '''
seqtk cutN -n 0 {input.fasta[0]} |\
calN50.js -L 3G - > {output}
        '''

rule calculate_gene_completeness:
    input:
        fasta = rules.panSN_renaming.output['fasta']
    output:
        metrics = expand('analyses/QC/completeness/compleasm_{{sample}}/{result}',result=('summary.txt','full_table.tsv'))
    params:
        _dir = lambda wildcards, output: PurePath(output['metrics'][0]).parent
    conda: 'compleasm'
    threads: 4
    resources:
        mem_mb_per_cpu = 12500,
        runtime = '2h'
    shell:
        '''
compleasm run -a {input.fasta[0]} -o {params._dir} -l cetartiodactyla -t {threads}

cut -f 3,2,1,11,9 {params._dir}/cetartiodactyla_odb10/full_table.tsv | sed '1d' > {output.metrics[1]}

rm -rf {params._dir}/cetartiodactyla_odb10
        '''

rule minimap2_reference_aligned:
    input:
        fasta = rules.panSN_renaming.output['fasta'],
        reference = config['reference']
    output:
        'analyses/QC/reference_alignment/{sample}.ARS_UCD2.0.paf.gz'
    threads: 4
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '2h'
    shell:
        '''
minimap2 -t {threads} --cs -cxasm10 {input.reference} {input.fasta[0]} |\
pigz -p {threads} -c > {output}
        '''

rule calculate_variant_level:
    input:
        paf = rules.minimap2_reference_aligned.output,
        reference = config['reference']
    output:
        vcf = multiext('analyses/QC/reference_alignment/{sample}.vcf.gz','','.csi')
    threads: 1
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '1h'
    shell:
        '''
pigz -p 2 -dc {input.paf} |\
sort -k6,6 -k8,8n |\
paftools.js call -f {input.reference} -s {wildcards.sample} - |\
bcftools view --write-index -o {output.vcf[0]}
        '''

rule calculate_reference_coverage:
    input:
        paf = rules.minimap2_reference_aligned.output,
        fai = config['reference'] + ".fai"
    output:
        bed = 'analyses/QC/reference_alignment/{sample}.covered.bed'
    threads: 1
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '30m'
    shell:
        '''
pigz -p 2 -dc {input.paf} |\
cut -f 6,8,9 |\
bedtools sort -faidx {input.fai} -i /dev/stdin |\
bedtools merge -d 0 -i /dev/stdin |\
bedtools genomecov -g {input.fai} -i /dev/stdin |\
awk -v OFS='\\t' '$2==0 {{A[$1]; U[$1]=$5; next}} {{A[$1];C[$1]=$5}} END {{for (k in C) {{print k,C[k]?C[k]:0,U[k]?U[k]:0}} }}' |\
sort -k1,1V > {output.bed}
        '''

rule summarise_sample_metrics:
    input:
        N50 = rules.calculate_N50.output,
        completeness = rules.calculate_gene_completeness.output['metrics'],
        variants = rules.calculate_variant_level.output['vcf'],
        bed = rules.calculate_reference_coverage.output['bed'],
        busco_map = config['busco_map']
    output:
        csv = 'analyses/QC/summary/{sample}.csv',
        busco = 'analyses/QC/completeness/{sample}.csv'
    resources:
        runtime = '10m'
    shell:
        '''
echo -n "{wildcards.sample}," > {output.csv}

awk '$1~/(SZ|NN)/ {{printf $2",";next}} {{if ($1=="NL"&&$2==50) {{printf $3","}} }}' {input.N50} >> {output.csv}

awk 'NR==FNR {{loc[$1]=$2;next}} {{++C[loc[$1]][$2]}} END {{ for (k in C) {{ print k,C[k]["Single"]?C[k]["Single"]:0,C[k]["Duplicated"]?C[k]["Duplicated"]:0,C[k]["Missing"]?C[k]["Missing"]:0 }} }}' {input.busco_map} {input.completeness[1]} > {output.busco}

awk '$1~/[[:digit:]]/ {{s+=$2;d+=$3;m+=$4;next}} {{S[$1]=$2;D[$1]=$3;M[$1]=$4}} END {{printf s","d","m","S["X"]","D["X"]","M["X"]","S["Y"]","D["Y"]","M["Y"]","}}' {output.busco} >> {output.csv}

bcftools stats {input.variants[0]} | awk '$1=="SN"&&$5~/(SNPs|indels)/ {{printf $6","}}' >> {output.csv}

awk -v OFS=',' '$1~/^[[:digit:]]/ {{A+=$2;++n;next}} {{B[$1]=$2}} END {{print A/n,B["X"]?B["X"]:0,B["Y"]?B["Y"]:0,B["MT"]?B["MT"]:0}}' {input.bed} >> {output.csv}
        '''

#bcftools merge --threads 2 -W -o {output.vcf[0]} {input.vcf}
#bcftools stats -r $(echo {1..29} | tr ' ' ',') -s - 67_samples.vcf.gz | grep "PSC"
rule summarise_all_metrics:
    input:
        metrics = lambda wildcards: expand(rules.summarise_sample_metrics.output['csv'],sample=determine_pangenome_samples(wildcards.graph)),
        vcf = lambda wildcards: expand(rules.calculate_variant_level.output['vcf'][0],sample=determine_pangenome_samples(wildcards.graph))
    output:
        metrics = 'analyses/QC_summary.{graph}.csv',
        #vcf = multiext('analyses/QC_variants.vcf.gz','','.csi')
    localrule: True
    shell:
        '''
{{ echo "sample,genome size,N contigs,NG50,autosome single copy,autosome duplicated copy,autosome missing copy,X single copy,X duplicated copy,X missing copy,Y single copy,Y duplicated copy,Y missing copy,SNPs,InDels,autosomes covered,X covered,Y covered,MT covered" ;  cat {input.metrics} ; }} > {output.metrics}
        '''
