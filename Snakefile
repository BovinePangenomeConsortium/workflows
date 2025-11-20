import polars as pl
import re


ALL_CHROMOSOME = list(map(str, range(1, 30))) + ["X", "Y", "MT"]

metadata = pl.read_csv(config["metadata"], infer_schema_length=10000).with_columns(
    pl.concat_str([pl.col("Animal ID"), pl.col("Haplotype")], separator="_").alias(
        "Filename"
    )
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
    graph="|".join(["every"] + graph_choices),
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
        match wildcards.graph:
            case (  # can I capture all cases based on graph_choices?
                "Cattle_n1" | "Cattle" | "All"
            ):  # Note does not include QC fails in "all"
                subset = metadata.filter(
                    pl.col("Graphs").str.split(";").list.contains(wildcards.graph)
                )
            case "breed_representative":
                subset = metadata.filter(
                    (pl.col("Species") == "Bos taurus")
                    & (pl.col("Breed representative") == "Y")
                )
            case "subspecies_representative":
                subset = metadata.filter(
                    (pl.col("Species") == "Bos taurus")
                    & (pl.col("Subspecies representative") == "Y")
                )
            case "every":
                subset = metadata
            case _:
                subset = metadata
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


# include: 'snakepit/pggb.smk'

# Pangenome graph analyses
# include: 'snakepit/pangenome_alignment.smk'


rule all:
    input:
        "analyses/QC_summary.every.csv",
        "analyses/satellites/every.txt",
        "analyses/pangene/every.ARS_UCD2.0_Ensembl.clustered.autosomes.tsv",
        "analyses/minigraph/every/L50/mg.gfa",
        "cactus/every/genome.full.unchopped.gfa.gz",
