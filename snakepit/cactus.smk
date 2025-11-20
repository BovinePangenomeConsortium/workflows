# very stripped back version of https://github.com/harvardinformatics/cactus-snakemake/blob/main/cactus_minigraph.smk
# TODO: optimize resource usage and add more config options


rule cactus_make_seqfile:
    input:
        metadata=config["metadata"],
    output:
        seqfile="cactus/{graph}/seqfile.tsv",
    localrule: True
    run:
        metadata = pl.read_csv(input.metadata, infer_schema_length=10000)
        if wildcards.graph == "every":
            subset = metadata
        else:
            subset = metadata.filter(
                pl.col("Graphs").str.split(";").list.contains(wildcards.graph)
            )
        with open(output.seqfile, "w") as fout:
            for row in subset.iter_rows(named=True):
                fout.write(
                    f"{row["Animal ID"]}.{row["Haplotype"]}\tdata/agc_assemblies/{row["Animal ID"]}_{row["Haplotype"]}.fa.gz\n"
                )


rule cactus_minigraph:
    input:
        seqfile=ancient(rules.cactus_make_seqfile.output["seqfile"]),
        fasta=expand(
            rules.panSN_renaming.output["fasta"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        gfa="cactus/{graph}/sv.gfa",
    params:
        reference_ID=config.get("reference_ID", "ARS_UCD2.0"),
    threads: 8
    resources:
        mem_mb_per_cpu=10000,
        runtime="4h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-minigraph \
$TMPDIR/minigraph \
{input.seqfile} \
{output.gfa} \
--reference {params.reference_ID}
"""


# Default cactus params use 6 threads per minigraph job
rule cactus_graphmap:
    input:
        seqfile=ancient(rules.cactus_make_seqfile.output["seqfile"]),
        gfa=rules.cactus_minigraph.output.gfa,
    output:
        multiext(
            "cactus/{graph}/graphmap",
            paf=".paf",
            fasta=".sv.gfa.fa",
            gaf=".gaf.gz",
            paf_filter_log=".paf.filter.log",
            paf_unfiltered=".paf.unfiltered.gz",
        ),
    params:
        reference_ID=config.get("reference_ID", "ARS_UCD2.0"),
    threads: 18
    resources:
        mem_mb_per_cpu=5000,
        runtime="4h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-graphmap \
$TMPDIR/graphmap \
{input.seqfile} \
{input.gfa} \
{output.paf} \
--outputFasta {output.fasta} \
--reference {params.reference_ID}
"""


checkpoint cactus_split:
    input:
        seqfile=ancient(rules.cactus_make_seqfile.output["seqfile"]),
        gfa=rules.cactus_minigraph.output.gfa,
        paf=rules.cactus_graphmap.output.paf,
    output:
        split_dir=directory("cactus/{graph}/split"),
        chromfile="cactus/{graph}/split/chromfile.txt",
    params:
        reference_ID=config.get("reference_ID", "ARS_UCD2.0"),
    threads: 2
    resources:
        mem_mb_per_cpu=20000,
        runtime="4h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-graphmap-split \
$TMPDIR/split \
{input.seqfile} \
{input.gfa} \
{input.paf} \
--outDir {output.split_dir} \
--reference {params.reference_ID}
"""


rule cactus_align:
    input:
        seqfile="cactus/{graph}/split/seqfiles/{contig}.seqfile",
        paf="cactus/{graph}/split/{contig}/{contig}.paf",
    output:
        multiext("cactus/{graph}/align/{contig}", hal=".hal", vg=".vg"),
    params:
        reference_ID=config.get("reference_ID", "ARS_UCD2.0"),
    threads: 4
    resources:
        mem_mb_per_cpu=10000,
        runtime="4h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-align \
$TMPDIR/align \
{input.seqfile} \
{input.paf} \
{output.hal} \
--pangenome \
--reference {params.reference_ID} \
--outVG
"""


def gather_alignments(wildcards):
    with checkpoints.cactus_split.get(**wildcards).output.chromfile.open() as f:
        return sorted(
            [f"cactus/{{graph}}/align/{row.split('\t')[0]}.vg" for row in f],
            key=lambda s: [
                int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", s)
            ],
        )


# TODO: add VCF outputs
rule cactus_join:
    input:
        vg=gather_alignments,
    output:
        multiext(
            "cactus/{graph}/genome",
            gfa_full=".full.unchopped.gfa.gz",
            gfa_clip=".clip.unchopped.gfa.gz",
            gbz_full=".full.gbz",
            gbz_clip=".clip.gbz",
            stats=".stats.tgz",
        ),
    params:
        out_dir=lambda wildcards, output: Path(output[0]).parent,
        reference_ID=config.get("reference_ID", "ARS_UCD2.0"),
        references=" ".join(config.get("additional_references", [])),
    threads: 8
    resources:
        mem_mb_per_cpu=10000,
        runtime="4h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-graphmap-join \
$TMPDIR/join \
--vg {input.vg} \
--outDir {params.out_dir} \
--outName {wildcards.graph} \
--reference {params.reference_ID} {params.references} \
--unchopped-gfa full clip \
--gbz full clip \
--vcfwave \
--vcf full clip
"""
