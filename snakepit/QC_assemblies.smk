rule calculate_N50:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
    output:
        "analyses/QC/contiguity/{sample}.N50",
    params:
        min_gap=config.get("cut_gap_min_len", 10),
    resources:
        runtime="30m",
    shell:
        """
seqtk cutN -n {params.min_gap} {input.fasta[0]} |\
calN50.js -L 3G - > {output}
        """


rule compleasm_download:
    output:
        directory("analyses/QC/completeness/{lineage}_odb12"),
    localrule: True
    conda:
        "compleasm"
    shell:
        """
compleasm download {wildcards.lineage} -L {output}
        """


rule calculate_gene_completeness:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
        db=expand(rules.compleasm_download.output, lineage="artiodactyla"),
    output:
        metrics=expand(
            "analyses/QC/completeness/compleasm_{{sample}}/{result}",
            result=("summary.txt", "full_table.tsv"),
        ),
    params:
        _dir=lambda wildcards, output: Path(output["metrics"][0]).parent,
    conda:
        "compleasm"
    threads: 4
    resources:
        mem_mb_per_cpu=12500,
        runtime="4h",
    shell:
        """
compleasm run -a {input.fasta[0]} -o {params._dir} -l artiodactyla -L {input.db} -t {threads}

cut -f 3,2,1,11,9 {params._dir}/artiodactyla_odb12/full_table.tsv | sed '1d' > {output.metrics[1]}

rm -rf {params._dir}/artiodactyla_odb12
        """


rule align_assembly_to_reference:
    input:
        fasta=rules.panSN_renaming.output["fasta"][0],
        reference=lambda wildcards: config["references"][wildcards.reference],
        reference_fai=lambda wildcards: f"{config['references'][wildcards.reference]}.fai",
    output:
        paf="analyses/QC/reference_alignment/{sample}.{reference}.asm.paf",
    params:
        mm2_opt=config.get("qc_mm2_opt", "-cx asm10 --cs"),
    threads: 6
    resources:
        mem_mb_per_cpu=8000,
        runtime="4h",
    shell:
        """
minimap2 -t {threads} {params.mm2_opt} {input.reference} {input.fasta} > {output.paf}
        """


rule calculate_variant_level:
    input:
        paf=rules.align_assembly_to_reference.output["paf"],
        reference=lambda wildcards: config["references"][wildcards.reference],
    output:
        vcf=multiext(
            "analyses/QC/reference_alignment/{sample}.{reference}.vcf.gz", "", ".csi"
        ),
    threads: 1
    resources:
        mem_mb_per_cpu=7500,
        runtime="1h",
    shell:
        """
sort -k6,6 -k8,8n {input.paf} |\
paftools.js call -f {input.reference} -s {wildcards.sample} - |\
bcftools view --write-index -o {output.vcf[0]}
        """


rule bedtools_makewindows:
    input:
        fai=lambda wildcards: f"{config['references'][wildcards.reference]}.fai",
    output:
        bed="analyses/QC/variant_density/{reference}.{size}.bed",
    localrule: True
    shell:
        """
bedtools makewindows -g {input.fai} -w {wildcards.size} |\
sort -k1,1 -k2,2n > {output.bed}
        """


rule bedtools_coverage_variant:
    input:
        bed=rules.bedtools_makewindows.output["bed"],
        vcf=rules.calculate_variant_level.output["vcf"],
    output:
        csv="analyses/QC/variant_density/{sample}.{reference}.{size}.csv",
    threads: 1
    resources:
        mem_mb_per_cpu=7500,
        runtime="4h",
    shell:
        """
bcftools query -i 'type="snps"' -f '%CHROM\\t%POS' {input.vcf[0]} |\
awk -v OFS='\\t' '{{print $1,$2,$2+1}}' |\
sort -k1,1 -k2,2n |\
bedtools coverage -a {input.bed} -b /dev/stdin -counts -sorted |\
awk 'BEGIN {{print "{wildcards.sample}"}} {{print $4}}' > {output.csv}
        """


rule gather_variant_coverages:
    input:
        bed=rules.bedtools_makewindows.output["bed"],
        csv=expand(
            rules.bedtools_coverage_variant.output,
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        csv="analyses/QC/variant_density/{reference}.{size}.csv.gz",
    threads: 4
    resources:
        mem_mb_per_cpu=5000,
        runtime="24h",
    shell:
        """
paste <(awk -v OFS='\\t' 'BEGIN {{print "chromosome","start"}} {{print $1,$2}}' {input.bed}) {input.csv} | pigz -p {threads} -c > {output.csv}
        """


rule calculate_reference_coverage:
    input:
        paf=rules.align_assembly_to_reference.output["paf"],
        fai=lambda wildcards: f"{config['references'][wildcards.reference]}.fai",
    output:
        bed="analyses/QC/reference_alignment/{sample}.{reference}.covered.bed",
    threads: 1
    resources:
        mem_mb_per_cpu=2500,
        runtime="30m",
    shell:
        """
cut -f 6,8,9 {input.paf} |\
bedtools sort -faidx {input.fai} -i /dev/stdin |\
bedtools merge -d 0 -i /dev/stdin |\
bedtools genomecov -g {input.fai} -i /dev/stdin |\
awk -v OFS='\\t' '$2==0 {{A[$1]; U[$1]=$5; next}} {{A[$1];C[$1]=$5}} END {{for (k in C) {{print k,C[k]?C[k]:0,U[k]?U[k]:0}} }}' |\
sort -k1,1V > {output.bed}
        """


rule summarise_sample_metrics:
    input:
        N50=rules.calculate_N50.output,
        busco=rules.calculate_gene_completeness.output["metrics"][0],
        variants=expand(
            rules.calculate_variant_level.output["vcf"],
            reference=REFERENCE_ID,
            allow_missing=True,
        ),
        bed=expand(
            rules.calculate_reference_coverage.output["bed"],
            reference=REFERENCE_ID,
            allow_missing=True,
        ),
    output:
        csv="analyses/QC/summary/{sample}.csv",
    resources:
        runtime="20m",
    shell:
        """
echo -n "{wildcards.sample}," > {output.csv}

awk '$1~/(SZ|NN)/ {{printf $2",";next}} {{if ($1=="NL"&&$2==50) {{printf $3","}} }}' {input.N50} >> {output.csv}

awk 'NR>1&&NR<8 {{printf $2";"}}' {input.busco} | sed 's/;;$/,/' >> {output.csv}

bcftools stats {input.variants[0]} | awk '$1=="SN"&&$5~/(SNPs|indels)/ {{printf $6","}}' >> {output.csv}

awk -v OFS=',' '$1~/^[[:digit:]]/ {{A+=$2;++n;next}} {{B[$1]=$2}} END {{print A/n,B["X"]?B["X"]:0,B["Y"]?B["Y"]:0,B["MT"]?B["MT"]:0}}' {input.bed} >> {output.csv}
        """


rule summarise_all_metrics:
    input:
        metrics=expand(
            rules.summarise_sample_metrics.output["csv"],
            sample=determine_pangenome_samples,
        ),
    output:
        metrics="analyses/QC_summary.{graph}.csv",
    localrule: True
    shell:
        """
{{ echo "sample,genome size,N contigs,NG50,BUSCOs,SNPs,InDels,autosomes covered,X covered,Y covered,MT covered" ;  cat {input.metrics} ; }} > {output.metrics}
        """


rule merge_metadata_and_qc:
    input:
        manifest=config["metadata"],
        qc=rules.summarise_all_metrics.output["metrics"],
    output:
        tsv="analyses/QC/BPC_metadata_with_QC.{graph}.tsv",
    params:
        script=workflow.source_path("../scripts/manifest.py"),
        thresholds=json.dumps(config.get("qc_thresholds", {})),
    localrule: True
    shell:
        """
python3 {params.script} merge-qc \
  --manifest {input.manifest} \
  --qc {input.qc} \
  --output {output.tsv} \
  --thresholds-json '{params.thresholds}'
        """


rule bcftools_merge:
    input:
        vcf=expand(
            rules.calculate_variant_level.output["vcf"][0],
            sample=determine_pangenome_samples,
            reference=REFERENCE_ID,
        ),
        fasta=config["references"][REFERENCE_ID],
    output:
        vcf=multiext("analyses/QC/VCF/{graph}.vcf.gz", "", ".csi"),
    params:
        regions=",".join(ALL_CHROMOSOME),
    threads: 4
    resources:
        mem_mb_per_cpu=12500,
        runtime="4h",
    shell:
        """
bcftools merge --threads {threads} -Ou -r {params.regions} {input.vcf} |\
bcftools +setGT --threads {threads} -Ou - -- -t . -n 0 |\
bcftools norm -f {input.fasta} -Ou |\
bcftools sort -T $TMPDIR -Ou |\
bcftools norm --threads {threads} -Ou -m+any |\
bcftools +fill-tags --threads {threads} -Ou - -- -t AF |\
bcftools annotate --threads {threads} --set-id '%VKX' -W -o {output.vcf[0]}
        """

rule bcftools_stats:
    input:
        vcf=rules.bcftools_merge.output["vcf"]
    output:
        stats="analyses/QC/VCF/{graph}.stats.csv",
    params:
        chroms=" ".join(ALL_CHROMOSOME),
    threads: 2
    resources:
        mem_mb_per_cpu=2500,
        runtime="2h"
    shell:
        """
echo "sample,contig,nNonRefHom,nTransitions,nTransversions,nIndels,nSingletons" > {output.stats}

for R in {params.chroms}
  do
  bcftools stats --threads {threads} -s - -r $R {input.vcf[0]} | grep "^PSC" | awk -v OFS=',' -v R=$R '{{print $3,R,$5,$7,$8,$9,$11}}' >> {output.stats}
done
        """


rule plink_PCA:
    input:
        vcf=rules.bcftools_merge.output["vcf"],
    output:
        multiext(
            "analyses/QC/PCA/{graph}.{variants}", ".prune.in", ".prune.out", ".eigenval", ".eigenvec"
        ),
    params:
        prefix=lambda wildcards, output: Path(output[2]).with_suffix(""),
        variants=lambda wildcards: '--snps-only --max-alleles 2' if wildcards.variants=='biSNPs' else ''
    threads: 4
    resources:
        mem_mb_per_cpu=10000,
    shell:
        """
plink2 --threads {threads} --cow --vcf {input.vcf[0]} --indep-pairwise 100kb 0.8 --make-pgen --out {params.prefix} {params.variants}

plink2 --threads {threads} --cow --pfile {params.prefix} --pca --maf 0.1 --out {params.prefix} --extract {output[0]}
        """

#needs forked simd-sketch in progress
rule simd_sketch_PCA:
    input:
        fasta = expand(rules.panSN_renaming.output["fasta"],
                sample=determine_pangenome_samples,
                allow_missing=True)
    output:
        matrix = "analyses/QC/PCA/{graph}.matrix"
    threads: 6
    resources:
         mem_mb_per_cpu=2500,
    shell:
        """
simd-sketch matrix --alg bottom -s 100000 -j {threads} {input.fasta} > {output.matrix}
        """

#TODO: compress fasta? compress out?
rule repeat_masker:
    input:
        fasta = rules.panSN_renaming.output["fasta"],
        repeats = config.get('repeat_library',[])
    output:
        multiext("analyses/masking/{sample}.fa",masked=".masked",out=".out",cat=".cat.gz",tbl=".tbl")
    params:
        _dir=lambda wildcards, output: Path(output[0]).parent,
    threads: 12
    resources:
        mem_mb_per_cpu = 3500,
        runtime = '48h'
    shell:
        """
RepeatMasker -xsmall -pa {threads} -e rmblast -lib {input.repeats} -q -no_is {input.fasta[0]} -dir {params._dir}
        """

rule gather_masking:
    input:
        out = expand(rules.repeat_masker.output["out"],sample=determine_pangenome_samples)
    output:
        "analyses/masking/masking_statistics.csv.gz"
    localrule: True
    shell:
        """
cat {input} > {output}
        """
