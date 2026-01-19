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
