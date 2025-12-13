import polars as pl
import re


ALL_CHROMOSOME = list(map(str, range(1, 30))) + ["X", "Y", "MT"]

metadata = (pl.read_csv(config["metadata"], infer_schema_length=10000)
        .with_columns(pl.concat_str([pl.col("Animal ID"), pl.col("Haplotype")], separator="_").alias("Filename"))
)

graph_choices = (
    metadata.select(
        pl.col("Graphs")
        .drop_nulls()
        .str.split(";")
        .list.explode()  # Flatten the list
        .unique()
    )
    .get_column("Graphs")
    .to_list()
)


wildcard_constraints:
    sample=r"[\w+\.\-_]+",
    graph="|".join(graph_choices),
    chromosome="|".join(ALL_CHROMOSOME),
    L=r"\d+",
    reference="|".join(config.get("references", [])),
    peptides="|".join(config.get("peptides", [])),


try:
    ANNOTATED_GENOMES = (
        metadata.filter(pl.col("Reference annotation") == "Y")
        .get_column("Filename")
        .to_list()
    )
except:
    ANNOTATED_GENOMES = []

ALL_CHROMOSOME = list(map(str, range(1, 30))) + ["X", "Y", "MT"]

alignment_metadata = (
    pl.read_csv(config["alignment_metadata"], infer_schema_length=10000)
    if "alignment_metadata" in config
    else pl.DataFrame()
)


def determine_pangenome_subset(wildcards):
    try:
        subset = metadata.filter(
            pl.col("Graphs").str.split(";").list.contains(wildcards.graph)
        )
    except AttributeError:
        subset = metadata
    return subset


def determine_pangenome_samples(wildcards):
    subset = determine_pangenome_subset(wildcards)
    return subset.get_column("Filename").to_list()


# Assembly level analyses
include: "snakepit/pangenome.smk"
include: "snakepit/QC_assemblies.smk"
include: "snakepit/centomeres.smk"
# Pangenome graph builders
include: "snakepit/minigraph.smk"
include: "snakepit/pangene.smk"
include: "snakepit/cactus.smk"

# Pangenome alignment
include: "snakepit/alignment.smk"

# include: 'snakepit/pggb.smk'

# Pangenome graph analyses
# include: 'snakepit/pangenome_alignment.smk'


rule all:
    input:
        "analyses/QC_summary.All.csv",
        "analyses/QC/PCA/Bovinae.biSNPs.eigenvec",
        "analyses/satellites/Bovinae.csv",
        "analyses/pangene/Bovinae.ARS_UCD2.0_Ensembl.clustered.tsv",
        "analyses/minigraph/Bovina/L50/mg.gfa",
        "cactus/BosTaurus_n1/graph.gbz",
