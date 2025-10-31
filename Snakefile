import polars as pl


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
    graph=r"|".join(graph_choices),
    chromosome=r"|".join(ALL_CHROMOSOME),
    L=r"\d+",
    reference="|".join(config.get("references")),
    peptides="|".join(config.get("peptides")),


ANNOTATED_GENOMES = (
    metadata.filter(pl.col("Reference annotation") == "Y")
    .get_column("Filename")
    .to_list()
)

ALL_CHROMOSOME = list(map(str, range(1, 30))) + ["X", "Y", "MT"]

alignment_metadata = (
    pl.read_csv(config["alignment_metadata"], infer_schema_length=10000)
    if "alignment_metadata" in config
    else pl.DataFrame()
)


# TODO: need to rewrite to handle AGC IDs
def determine_pangenome_samples(wildcards):
    try:
        match wildcards.graph:
            case (
                "n=1 Cattle" | "Cattle" | "All"
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
            case "all":
                subset = metadata
            case _:
                subset = metadata
    except AttributeError:
        subset = metadata
    return subset.get_column("Filename").to_list()


include: "snakepit/pangenome.smk"
include: "snakepit/QC_assemblies.smk"


# include: 'snakepit/minigraph.smk'
# include: 'snakepit/pggb.smk'
# include: 'snakepit/pangene.smk'
# include: 'snakepit/pangenome_alignment.smk'
# include: 'snakepit/centomeres.smk'


rule all:
    input:
        "analyses/QC_summary.All.csv",
        #expand('analyses/pggb/medium_test/p95_s5000/{chromosome}.k31.POAasm20.unchop.graph.png',chromosome=ALL_CHROMOSOME)
