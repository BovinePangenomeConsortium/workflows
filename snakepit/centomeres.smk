rule KMC_count:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
    output:
        kmers=multiext(
            "analyses/satellites/{sample}.kmc", ".kmc_pre", ".kmc_suf", ".txt"
        ),
    params:
        prefix=lambda wildcards, output: Path(output[0]).with_suffix(""),
    threads: 2
    resources:
        mem_mb_per_cpu=10000,
        runtime="1h",
    shell:
        """
kmc -k151 -t{threads} -ci20 -cs100000 -fm {input.fasta[0]} {params.prefix} $TMPDIR

kmc_dump {params.prefix} {output.kmers[2]}
        """


rule SRF:
    input:
        kmers=rules.KMC_count.output["kmers"][2],
    output:
        satellites="analyses/satellites/{sample}.srf.fa",
    resources:
        mem_mb_per_cpu=2500,
        runtime="30m",
    shell:
        """
        srf -p {wildcards.sample} -l 40 -c 100 {input.kmers} > {output.satellites}
        """


rule minimap2_srf:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
        satellites=rules.SRF.output["satellites"],
    output:
        paf="analyses/satellites/{sample}.srf-aln.paf",
    threads: 6
    resources:
        mem_mb_per_cpu=2500,
        runtime="24h",
    shell:
        """
minimap2 -t {threads} -c -N1000000 -f1000 -r100,100 {input.satellites} {input.fasta[0]} > {output.paf}
        """


rule srfutils:
    input:
        paf=rules.minimap2_srf.output["paf"],
    output:
        bed="analyses/satellites/{sample}.srf-aln.bed",
        abun="analyses/satellites/{sample}.srf-aln.abun",
    threads: 1
    resources:
        mem_mb_per_cpu=5000,
        runtime="1h",
    shell:
        """
        srfutils.js paf2bed {input.paf} > {output.bed}

        srfutils.js bed2abun -g 3G {output.bed} > {output.abun}
        """


rule srf_gather:
    input:
        paf=expand(
            rules.srfutils.output["abun"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        gfa="analyses/satellites/{graph}.txt",
    localrule: True
    shell:
        """
        touch {output}
        """

# TODO: we can use https://github.com/richarddurbin/rotate and cd-hit-est to realign+cluster sequences
# Will require some careful thought on picking initial frame sequence for each major sat type, as well as pre-filtering low abundance repeats
# e.g. for BTSAT of length 1400
# rotate -s AAGGGGTTCCCGAC sats.fa > rotated.fa
# cd-hit-est -i rotated.fa -o r.fa -c 0.99 -n 11 -G 0 -aS 0.9 -d 0 -T 1
