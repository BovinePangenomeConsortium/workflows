rule vg_gbwt:
    input:
        gbz="analyses/cactus/{graph}/{graph}.gbz",
    output:
        ri="analyses/cactus/{graph}/{graph}.ri",
    threads: 2
    resources:
        mem_mb_per_cpu=30000,
        runtime="4h",
    shell:
        """
vg gbwt \
--gbz-input {input.gbz} \
--temp-dir $TMPDIR \
--r-index {output.ri} \
--num-threads {threads}
        """


rule vg_index:
    input:
        gbz="analyses/cactus/{graph}/{graph}.gbz",
    output:
        dist="analyses/cactus/{graph}/{graph}.dist",
    threads: 4
    resources:
        mem_mb_per_cpu=50000,
        runtime="4h",
    shell:
        """
vg index \
--temp-dir $TMPDIR \
--no-nested-distance \
--dist-name {output.dist} \
--threads {threads} \
{input.gbz}
        """


rule vg_haplotype:
    input:
        gbz="analyses/cactus/{graph}/{graph}.gbz",
        ri=rules.vg_gbwt.output["ri"],
        dist=rules.vg_index.output["dist"],
    output:
        hapl="analyses/cactus/{graph}/{graph}.hapl",
    threads: 4
    resources:
        mem_mb_per_cpu=40000,
        runtime="4h",
    shell:
        """
vg haplotypes \
--threads {threads} \
--haplotype-output {output.hapl} \
--distance-index {input.dist} \
--r-index {input.ri} \
{input.gbz}
        """

rule deacon_build:
    input:
        fasta=rules.panSN_renaming.output["fasta"]
    output:
        index="analyses/deacon/{sample}.idx"
    threads: 4
    resources:
        mem_mb_per_cpu=5000,
        runtime="1h",
    shell:
        """
deacon index build -k 29 -w 11 --threads {threads} {input.fasta[0]} > {output.index}
        """

rule deacon_pangenome:
    input:
        reference=expand(rules.deacon_build.output["index"],sample=config.get("reference_ID", "ARS_UCD2.0").replace(".","_")),
        indexes=expand(rules.deacon_build.output["index"],sample=determine_pangenome_samples),
    output:
        multiext("analyses/deacon/{graph}",index=".union",nonreference=".diff")
    threads: 1
    resources:
        mem_mb_per_cpu=50000,
        runtime="4h",
    shell:
        """
deacon index union {input.indexes} > {output.index}
deacon index diff {output.index} {input.reference} > {output.nonreference}
        """

rule seqtk_QC:
    input:
        fastq="/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/tvanelst/SRA/Bos/PE9/{sample}_R{N}.final.fq.gz",
    output:
        fasta="analyses/deacon/{sample}.raw.R{N}.fasta.gz",
    threads: 1
    resources:
        mem_mb_per_cpu=5000,
        runtime="1h",
    shell:
        """
seqtk seq -A -L 150 {input.fastq} > {output.fasta}
        """

#seqtk seq -A -L 150 {input.fastq[0]} > $TMPDIR/R1.fa
#seqtk seq -A -L 150 {input.fastq[1]} > $TMPDIR/R2.fa
rule deacon_filter:
    input:
        fastq=lambda wildcards: expand(rules.deacon_filter.output["fasta"],reads="unmapped",allow_missing=True) if wildcards.reads != 'unmapped' else expand("/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/tvanelst/SRA/Bos/PE9/{sample}_R{N}.final.fq.gz",N=(1,2),allow_missing=True),
        index=lambda wildcards: expand(rules.deacon_pangenome.output["nonreference"],graph=wildcards.reads) if wildcards.reads != 'unmapped' else expand(rules.deacon_build.output["index"],sample=config.get("reference_ID", "ARS_UCD2.0").replace(".","_")),
    output:
        summary="analyses/deacon/{sample}.{reads}.json",
        fasta=expand("analyses/deacon/{sample}.{reads}.R{N}.fasta.gz",N=(1,2),allow_missing=True),
    threads: 4
    resources:
        mem_mb_per_cpu=12000,
        runtime="4h",
    shell:
        """
deacon filter --threads {threads} --summary {output.summary} --quiet --deplete \
--output {output.fasta[0]} --output2 {output.fasta[1]} \
{input.index} {input.fastq}
        """
