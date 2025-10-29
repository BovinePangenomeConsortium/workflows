# very stripped back version of https://github.com/harvardinformatics/cactus-snakemake/blob/main/cactus_minigraph.smk
# TODO: optimize resource usage and add more config options
# TODO: check on split and gather steps

rule all:
    input:
        "cactus/all/all.gfa.gz"

rule cactus_copy_seqfile:
    input:
        seqfile=config["seq_file"]
    output:
        seqfile="cactus/{graph}/seqfile.txt"
    localrule: True
    shell: '''
cp {input.seqfile} {output.seqfile}
'''


# bottleneck is sequential
rule cactus_minigraph:
    input:
        seqfile=config["seq_file"]
    output:
        gfa="cactus/{graph}/sv.gfa"
    threads: 8
    resources:
        mem_mb_per_cpu=10000,
        runtime="4h"
    container: config["cactus_image"]
    shell: '''
cactus-minigraph \
$TMPDIR/minigraph \
{input.seqfile} \
{output.gfa} \
--reference {config[reference_ID]}
'''


# bottleneck is parallel
rule cactus_graphmap:
    input:
        seqfile=ancient(rules.cactus_copy_seqfile.output["seqfile"]),
        gfa=rules.cactus_minigraph.output.gfa
    output:
        multiext(
            "cactus/{graph}/graphmap",
            paf=".paf",
            fasta=".sv.gfa.fa",
            gaf=".gaf.gz",
            paf_filter_log=".paf.filter.log",
            paf_unfiltered=".paf.unfiltered.gz",
        )
    threads: 8
    resources:
        mem_mb_per_cpu=10000,
        runtime="4h"
    container: config["cactus_image"]
    shell: '''
cactus-graphmap \
$TMPDIR/graphmap \
{input.seqfile} \
{input.gfa} \
{output.paf} \
--outputFasta {output.fasta} \
--reference {config[reference_ID]}
'''


checkpoint cactus_split:
    input:
        seqfile=ancient(rules.cactus_copy_seqfile.output["seqfile"]),
        gfa=rules.cactus_minigraph.output.gfa,
        paf=rules.cactus_graphmap.output.paf
    output:
        split_dir=directory("cactus/{graph}/split"),
        chromfile="cactus/{graph}/split/chromfile.txt"
    threads: 8
    resources:
        mem_mb_per_cpu=10000,
        runtime="4h"
    container: config["cactus_image"]
    shell: '''
cactus-graphmap-split \
$TMPDIR/split \
{input.seqfile} \
{input.gfa} \
{input.paf} \
--outDir {output.split_dir} \
--reference {config[reference_ID]}
'''


# per chromosome bottleneck
rule cactus_align:
    input:
        seqfile="cactus/{graph}/split/seqfiles/{contig}.seqfile",
        paf="cactus/{graph}/split/{contig}/{contig}.paf"
    output:
        multiext("cactus/{graph}/align/{contig}", hal=".hal", vg=".vg")
    threads: 8
    resources:
        mem_mb_per_cpu=10000,
        runtime="4h"
    container: config["cactus_image"]
    shell: '''
cactus-align \
$TMPDIR/align \
{input.seqfile} \
{input.paf} \
{output.hal} \
--pangenome \
--reference {config[reference_ID]} \
--outVG
'''


import re


def gather_alignments(wildcards):
    with checkpoints.cactus_split.get(**wildcards).output.chromfile.open() as f:
        return sorted(
            [f"cactus/{{graph}}/align/{row.split('\t')[0]}.vg" for row in f],
            key=lambda s: [int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", s)],
        )


# TODO: upgrade other references here?
rule cactus_join:
    input:
        vg=gather_alignments
    output:
        multiext(
            "cactus/{graph}/{graph}",
            hal=".full.hal",
            gfa=".gfa.gz",
            vcf=".vcf.gz",
            vcf_index=".vcf.gz.tbi",
            gbz=".gbz",
            raw_vcf=".raw.vcf.gz",
            raw_vcf_index=".raw.vcf.gz.tbi",
            stats=".stats.tgz",
        )
    params:
        out_dir=lambda wildcards, output: Path(output[0]).parent
    threads: 8
    resources:
        mem_mb_per_cpu=10000,
        runtime="4h"
    container: config["cactus_image"]
    shell: '''
cactus-graphmap-join \
$TMPDIR/join \
--vg {input.vg} \
--outDir {params.out_dir} \
--outName {wildcards.graph} \
--reference {config[reference_ID]} \
--unchopped-gfa full clip \
--gbz full clip \
--vcfwave \
--vcf full clip
'''