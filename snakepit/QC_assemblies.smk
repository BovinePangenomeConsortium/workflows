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

rule compleasm_download:
    output:
        directory('analyses/QC/completeness/{lineage}_odb12')
    localrule: True
    conda: 'compleasm'
    shell:
        '''
        compleasm download {wildcards.lineage} -L {output}
        '''

rule calculate_gene_completeness:
    input:
        fasta = rules.panSN_renaming.output['fasta'],
        db = expand(rules.compleasm_download.output,lineage='artiodactyla')
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
compleasm run -a {input.fasta[0]} -o {params._dir} -l artiodactyla -t {threads}

cut -f 3,2,1,11,9 {params._dir}/artiodactyla_odb12/full_table.tsv | sed '1d' > {output.metrics[1]}

rm -rf {params._dir}/artiodactyla_odb12
        '''

rule minimap2_reference_aligned:
    input:
        fasta = rules.panSN_renaming.output['fasta'],
        reference = lambda wildcards: config['references'][wildcards.reference]
    output:
        'analyses/QC/reference_alignment/{sample}.{reference}.paf.gz'
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
        paf = rules.ragtag_scaffold.output[5],
        #paf = rules.minimap2_reference_aligned.output,
        reference = lambda wildcards: config['references'][wildcards.reference]
    output:
        vcf = multiext('analyses/QC/reference_alignment/{sample}.{reference}.vcf.gz','','.csi')
    threads: 1
    resources:
        mem_mb_per_cpu = 7500,
        runtime = '1h'
    shell:
        '''
sort -k6,6 -k8,8n {input.paf} |\
paftools.js call -f {input.reference} -s {wildcards.sample} - |\
bcftools view --write-index -o {output.vcf[0]}
        '''

rule bedtools_makewindows:
    input:
        fai = lambda wildcards: f"{config['references'][wildcards.reference]}.fai"
    output:
        bed = 'analyses/QC/variant_density/{reference}.{size}.bed'
    localrule: True
    shell:
        '''
bedtools makewindows -g {input.fai} -w {wildcards.size} |\
sort -k1,1 -k2,2n > {output.bed}
        '''

rule bedtools_coverage_variant:
    input:
        bed = rules.bedtools_makewindows.output['bed'],
        vcf = rules.calculate_variant_level.output['vcf']
    output:
        csv = 'analyses/QC/variant_density/{sample}.{reference}.{size}.csv'
    threads: 1
    resources:
        mem_mb_per_cpu = 7500,
        runtime = '4h'
    shell:
        '''
bcftools query -i 'type="snps"' -f '%CHROM\\t%POS' {input.vcf[0]} |\
awk -v OFS='\\t' '{{print $1,$2,$2+1}}' |\
sort -k1,1 -k2,2n |\
bedtools coverage -a {input.bed} -b /dev/stdin -counts -sorted |\
awk 'BEGIN {{print "{wildcards.sample}"}} {{print $4}}' > {output.csv}
        '''

rule gather_variant_coverages:
    input:
        bed = rules.bedtools_makewindows.output['bed'],
        csv = expand(rules.bedtools_coverage_variant.output,sample=determine_pangenome_samples,allow_missing=True)
    output:
        csv = 'analyses/QC/variant_density/{reference}.{size}.csv.gz'
    threads: 4
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '24h'
    shell:
        '''
paste <(awk -v OFS='\\t' 'BEGIN {{print "chromosome","start"}} {{print $1,$2}}' {input.bed}) {input.csv} | pigz -p {threads} -c > {output.csv}
        '''

rule calculate_reference_coverage:
    input:
        paf = rules.ragtag_scaffold.output[5],
        #paf = rules.minimap2_reference_aligned.output,
        fai = lambda wildcards: f"{config['references'][wildcards.reference]}.fai"
    output:
        bed = 'analyses/QC/reference_alignment/{sample}.{reference}.covered.bed'
    threads: 1
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '30m'
    shell:
        '''
cut -f 6,8,9 {input.paf} |\
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
        variants = expand(rules.calculate_variant_level.output['vcf'],reference='ARS_UCD2.0',allow_missing=True),
        bed = expand(rules.calculate_reference_coverage.output['bed'],reference='ARS_UCD2.0',allow_missing=True),
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
        metrics = expand(rules.summarise_sample_metrics.output['csv'],sample=determine_pangenome_samples),
        #vcf = expand(rules.calculate_variant_level.output['vcf'][0],sample=determine_pangenome_samples,reference='ARS-UCD2.0')
    output:
        metrics = 'analyses/QC_summary.{graph}.csv',
        #vcf = multiext('analyses/QC_variants.vcf.gz','','.csi')
    localrule: True
    shell:
        '''
{{ echo "sample,genome size,N contigs,NG50,autosome single copy,autosome duplicated copy,autosome missing copy,X single copy,X duplicated copy,X missing copy,Y single copy,Y duplicated copy,Y missing copy,SNPs,InDels,autosomes covered,X covered,Y covered,MT covered" ;  cat {input.metrics} ; }} > {output.metrics}
        '''

rule bcftools_merge:
    input:
        vcf = expand(rules.calculate_variant_level.output['vcf'][0],sample=determine_pangenome_samples,reference='ARS_UCD2.0')
    output:
        vcf = multiext('analyses/{graph}.vcf.gz','','.csi'),
        stats = 'analyses/{graph}.vcf.stats'
    threads: 4
    resources:
        mem_mb_per_cpu = 12500,
        walltime = '4h'
    shell:
        '''
        bcftools merge --threads {threads} -Ou -r $(echo {{1..29}} | tr ' ' ',') {input.vcf} |\
        bcftools +setGT --threads {threads} -Ou - -- -t . -n 0 |\
        bcftools norm --threads {threads} -Ou -m-any |\
        bcftools +fill-tags --threads {threads} -Ou - -- -t AF |\
        bcftools annotate --threads {threads} --set-id '%VKX' -W -o {output.vcf[0]}

        bcftools stats --threads {threads} -s - {output.vcf[0]} | grep "^PSC" | cut -f 3,5,7-9,11 > {output.stats}
        '''

rule plink_PCA:
    input:
        vcf = rules.bcftools_merge.output['vcf']
    output:
        multiext('analyses/PCA/{graph}','.prune.in','.prune.out','.eigenval','.eigenvec')
    params:
        prefix = lambda wildcards, output: Path(output[2]).with_suffix('')
    threads: 4
    resources:
        mem_mb_per_cpu = 10000
    shell:
        '''
        plink2 --threads {threads} --cow --vcf {input.vcf[0]} --indep-pairwise 100kb 0.8 --make-pgen --out {params.prefix} --snps-only --max-alleles 2

        plink2 --threads {threads} --cow --pfile {input.vcf[0]} --pca --maf 0.1 --out {params.prefix} --extract {output[0]}
        '''
