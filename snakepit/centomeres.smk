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


# Note: we set `-l 40` since the smallest expected bovine SAT is 46 bp, and to avoid picking up the 23 bp sub-motifs in other SATs
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
        runtime="4h",
    shell:
        """
minimap2 -t {threads} -c -N1000000 -f1000 -r100,100 <(srfutils.js enlong {input.satellites}) {input.fasta[0]} > {output.paf}
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
        abun=expand(
            rules.srfutils.output["abun"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
        satellites=expand(
            rules.SRF.output["satellites"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        csv="analyses/satellites/{graph}.csv",
    localrule: True
    run:
        srf_pattern = re.compile(r'>(?P<sample>\w+)#circ(?P<contig>\d+)-(?P<length>\d+) min=(?P<min>\d+),max=(?P<max>\d+),avg=(?P<avg>\d+)')
        rows = []
        for fname in input.satellites:
            with open(fname) as fin:
                for line in fin:
                    if line[0] == '>':
                        rows.append(srf_pattern.match(line.rstrip()).groupdict())
                    else:
                        rows[-1]['fasta']=line.rstrip()
        satellites = pl.DataFrame(rows)

        abun_pattern = re.compile(r'(?P<sample>\w+)#circ(?P<contig>\d+)-(?P<length>\d+)\t(?P<bp>\d+)\t(?P<mean_identity>(NaN|\d*\.?\d+))\t(?P<filtered>\d*\.?\d+)\t(?P<unfiltered>\d*\.?\d+)')

        rows = []
        for fname in input.abun:
            with open(fname) as fin:
                for line in fin:
                    rows.append(abun_pattern.match(line.rstrip()).groupdict())

        abundances = pl.DataFrame(rows)
        (satellites.join(abundances,on=['sample','contig','length'])
                   .with_columns(pl.col('mean_identity').str.replace("NaN","-1"))
                   .write_csv(output['csv'])
        )


# TODO: we can use https://github.com/richarddurbin/rotate and cd-hit-est to realign+cluster sequences
# Will require some careful thought on picking initial frame sequence for each major sat type, as well as pre-filtering low abundance repeats
# e.g. for BTSAT of length 1400
# rotate -s AAGGGGTTCCCGAC sats.fa > rotated.fa
# cd-hit-est -i rotated.fa -o r.fa -c 0.99 -n 11 -G 0 -aS 0.9 -d 0 -T 1
