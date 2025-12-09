# very stripped back version of https://github.com/harvardinformatics/cactus-snakemake/blob/main/cactus_minigraph.smk

rule cactus_make_seqfile:
    input:
        fasta=expand(rules.panSN_renaming.output["fasta"],sample=determine_pangenome_samples,allow_missing=True),
        metadata=config["metadata"],
    output:
        seqfile="analyses/cactus/{graph}/seqfile.tsv",
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
                    f"{row["Animal ID"]}.{row["Haplotype"]}\t{Path(config["assemblies_dir"]) / f"{row["Animal ID"]}_{row["Haplotype"]}.fa.gz"}\n"
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
        gfa="analyses/cactus/{graph}/sv.gfa",
    params:
        reference_ID=config.get("reference_ID", "ARS_UCD2.0"),
        jobstore=lambda wildcards,output: (Path("$TMPDIR") if config.get('jobstore',True) else Path(output.gfa).parent) / "minigraph_TEMPDIR",
    threads: 18
    resources:
        mem_mb_per_cpu=12000,
        runtime="120h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-minigraph \
{params.jobstore} \
--workDir {params.jobstore} \
{input.seqfile} \
{output.gfa} \
--reference {params.reference_ID}
"""


# Default cactus params use 6 threads per minigraph job, change with --mapCores
rule cactus_graphmap:
    input:
        seqfile=ancient(rules.cactus_make_seqfile.output["seqfile"]),
        gfa=rules.cactus_minigraph.output.gfa,
    output:
        multiext(
            "analyses/cactus/{graph}/graphmap",
            paf=".paf",
            fasta=".sv.gfa.fa",
            gaf=".gaf.gz",
            paf_filter_log=".paf.filter.log",
            paf_unfiltered=".paf.unfiltered.gz",
        ),
    params:
        reference_ID=config.get("reference_ID", "ARS_UCD2.0"),
        jobstore=lambda wildcards,output: (Path("$TMPDIR") if config.get('jobstore',True) else Path(output.paf).parent) / "graphmap_TEMPDIR",
    threads: 18
    resources:
        mem_mb_per_cpu=9000,
        runtime="120h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-graphmap \
{params.jobstore} \
--workDir {params.jobstore} \
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
        split_dir=directory("analyses/cactus/{graph}/split"),
        chromfile="analyses/cactus/{graph}/split/chromfile.txt",
    params:
        reference_ID=config.get("reference_ID", "ARS_UCD2.0"),
        jobstore=lambda wildcards,output: (Path("$TMPDIR") if config.get('jobstore',True) else Path(output.chromfile).parent) / "split_TEMPDIR",
    threads: 4
    resources:
        mem_mb_per_cpu=20000,
        runtime="48h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-graphmap-split \
{params.jobstore} \
--workDir {params.jobstore} \
{input.seqfile} \
{input.gfa} \
{input.paf} \
--outDir {output.split_dir} \
--reference {params.reference_ID}
"""


rule cactus_align:
    input:
        seqfile="analyses/cactus/{graph}/split/seqfiles/{contig}.seqfile",
        paf="analyses/cactus/{graph}/split/{contig}/{contig}.paf",
    output:
        multiext("analyses/cactus/{graph}/align/{contig}", hal=".hal", vg=".vg"),
    params:
        reference_ID=config.get("reference_ID", "ARS_UCD2.0"),
        jobstore=lambda wildcards,output: (Path("$TMPDIR") if config.get('jobstore',True) else Path(output.vg).parent) / "align_TEMPDIR",
    threads: 6
    resources:
        mem_mb_per_cpu=15000,
        runtime="48h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-align \
{params.jobstore} \
--workDir {params.jobstore} \
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
            [f"analyses/cactus/{{graph}}/align/{row.split('\t')[0]}.vg" for row in f],
            key=lambda s: [
                int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", s)
            ],
        )


# TODO: add VCF outputs
# split out long jobs
rule cactus_join:
    input:
        vg=gather_alignments,
    output:
        multiext(
            "analyses/cactus/{graph}/{graph}",
            gfa_full=".full.unchopped.gfa.gz",
            gbz_full=".full.gbz",
            snarls_full=".full.snarls",
            gfa_clip=".unchopped.gfa.gz",
            gbz_clip=".gbz",
            snarls_clip=".snarls",
            stats=".stats.tgz",
        ),
    params:
        out_dir=lambda wildcards, output: Path(output[0]).parent,
        reference_ID=config.get("reference_ID", "ARS_UCD2.0"),
        references=" ".join(config.get("additional_references", [])),
        jobstore=lambda wildcards,output: (Path("$TMPDIR") if config.get('jobstore',True) else Path(output.stats).parent) / "join_TEMPDIR",
    threads: 16
    resources:
        mem_mb_per_cpu=10000,
        runtime="24h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-graphmap-join \
{params.jobstore} \
--workDir {params.jobstore} \
--vg {input.vg} \
--outDir {params.out_dir} \
--outName {wildcards.graph} \
--reference {params.reference_ID} {params.references} \
--unchopped-gfa full clip \
--gbz full clip 
"""
#--vcf full clip \
#--vcfwave

