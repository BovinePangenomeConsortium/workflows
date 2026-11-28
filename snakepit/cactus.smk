# very stripped back version of https://github.com/harvardinformatics/cactus-snakemake/blob/main/cactus_minigraph.smk

rule cactus_prepare_fasta:
    input:
        fasta=rules.panSN_renaming.output["fasta"][0],
    output:
        fasta=temp(
            multiext("analyses/cactus/inputs/{sample}.fa.gz", "", ".fai", ".gzi")
        ),
    params:
        script=workflow.source_path("../scripts/format_pansn.py"),
    threads: 2
    resources:
        mem_mb_per_cpu=2000,
        runtime="30m",
    shell:
        """
python3 {params.script} strip-pansn \
  --fasta {input.fasta} \
  --output {output.fasta[0]} \
  --threads {threads}
        """


rule cactus_make_seqfile:
    input:
        fasta=expand(
            rules.cactus_prepare_fasta.output["fasta"][0],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
        metadata=config["metadata"],
    output:
        seqfile="analyses/cactus/{graph}/seqfile.tsv",
    localrule: True
    run:
        samples = determine_pangenome_samples(wildcards)
        ref_ids = [
            config.get("reference_ID", REFERENCE_ID),
            *config.get("additional_references", []),
        ]
        with open(output.seqfile, "w") as fout:
            for sample in samples:
                rec = manifest_records[sample]
                seq_name = _bpc_manifest.cactus_seqfile_name(rec, ref_ids)
                fout.write(
                    f"{seq_name}\tanalyses/cactus/inputs/{rec.sample}.fa.gz\n"
                )


rule cactus_minigraph:
    input:
        seqfile=ancient(rules.cactus_make_seqfile.output["seqfile"]),
        fasta=expand(
            rules.cactus_prepare_fasta.output["fasta"][0],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        gfa="analyses/cactus/{graph}/sv.gfa",
    params:
        reference_ID=config.get("reference_ID", REFERENCE_ID),
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
{input.seqfile} \
{output.gfa} \
--reference {params.reference_ID}
"""


# Default cactus params use 6 threads per minigraph job, change with --mapCores
rule cactus_graphmap:
    input:
        seqfile=ancient(rules.cactus_make_seqfile.output["seqfile"]),
        fasta=expand(
            rules.cactus_prepare_fasta.output["fasta"][0],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
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
        reference_ID=config.get("reference_ID", REFERENCE_ID),
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
{input.seqfile} \
{input.gfa} \
{output.paf} \
--outputFasta {output.fasta} \
--reference {params.reference_ID}
"""


checkpoint cactus_split:
    input:
        seqfile=ancient(rules.cactus_make_seqfile.output["seqfile"]),
        fasta=expand(
            rules.cactus_prepare_fasta.output["fasta"][0],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
        gfa=rules.cactus_minigraph.output.gfa,
        paf=rules.cactus_graphmap.output.paf,
    output:
        split_dir=directory("analyses/cactus/{graph}/split"),
        chromfile="analyses/cactus/{graph}/split/chromfile.txt",
    params:
        reference_ID=config.get("reference_ID", REFERENCE_ID),
        jobstore=lambda wildcards,output: (Path("$TMPDIR") if config.get('jobstore',True) else Path(output.chromfile).parent) / "split_TEMPDIR",
    threads: 4
    resources:
        mem_mb_per_cpu=40000,
        runtime="48h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-graphmap-split \
{params.jobstore} \
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
        reference_ID=config.get("reference_ID", REFERENCE_ID),
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
            "analyses/cactus/{graph}/graph",
            gfa_full=".full.unchopped.gfa.gz",
            #gbz_full=".full.gbz",
            #snarls_full=".full.snarls",
            gfa_clip=".unchopped.gfa.gz",
            gbz_clip=".gbz",
            hapl_clip=".hapl",
            snarls_clip=".snarls",
            stats=".stats.tgz",
        ),
    params:
        out_dir=lambda wildcards, output: Path(output[0]).parent,
        reference_ID=config.get("reference_ID", REFERENCE_ID),
        references=" ".join(config.get("additional_references", [])),
        jobstore=lambda wildcards,output: (Path("$TMPDIR") if config.get('jobstore',True) else Path(output.stats).parent) / "join_TEMPDIR",
        _workdir=lambda wildcards,output: (Path("$TMPDIR") if config.get('jobstore',True) else Path(output.stats).parent),
    threads: 12
    resources:
        mem_mb_per_cpu=30000,
        runtime="120h",
    container:
        config.get("cactus_image", "docker://docker.io/cactus/cactus:latest")
    shell:
        """
cactus-graphmap-join \
{params.jobstore} \
--workDir {params._workdir} \
--vg {input.vg} \
--outDir {params.out_dir} \
--outName graph \
--reference {params.reference_ID} {params.references} \
--unchopped-gfa full clip \
--snarlStats full clip \
--gbz clip \
--haplo clip \
--giraffe clip \
--lrGiraffe clip
"""
#--vcf full clip \
#--vcfwave
#https://github.com/Han-Cao/collapse-bubble
