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
        mem_mb_per_cpu=8000,
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
