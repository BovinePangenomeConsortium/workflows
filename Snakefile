import importlib.util
import json
from pathlib import Path
import polars as pl
import re

_manifest_spec = importlib.util.spec_from_file_location(
    "bpc_manifest", workflow.source_path("scripts/manifest.py")
)
_bpc_manifest = importlib.util.module_from_spec(_manifest_spec)
_manifest_spec.loader.exec_module(_bpc_manifest)

ALL_CHROMOSOME = [str(c) for c in config.get("chromosomes", _bpc_manifest.DEFAULT_CHROMOSOMES)]
REFERENCE_ID = config.get(
    "reference_ID",
    next(iter(config.get("references", {"ARS_UCD2.0": ""})), "ARS_UCD2.0"),
)
RAW_ASSEMBLIES_DIR = Path(config.get("raw_assemblies_dir", "data/raw_assemblies"))
STAGED_ASSEMBLIES_DIR = Path(config.get("staged_assemblies_dir", "data/staged_assemblies"))
ASSEMBLIES_DIR = Path(config.get("assemblies_dir", "data/assemblies"))

manifest_list = _bpc_manifest.load_manifest(config["metadata"], config=config)
manifest_records = {r.sample: r for r in manifest_list}
metadata = pl.DataFrame([r.to_metadata_row() for r in manifest_list])

graph_choices = _bpc_manifest.resolve_graph_choices(manifest_list, config)


wildcard_constraints:
    sample=r"[\w+\.\-_]+",
    graph="|".join(graph_choices),
    chromosome="|".join(ALL_CHROMOSOME),
    L=r"\d+",
    reference="|".join(config.get("references", [])),
    peptides="|".join(config.get("peptides", [])),


ANNOTATED_GENOMES = [
    r.sample
    for r in manifest_list
    if r.reference_annotation or r.format_mode == "ragtag_preserve"
]

alignment_metadata = (
    pl.read_csv(config["alignment_metadata"], infer_schema_length=10000)
    if "alignment_metadata" in config
    else pl.DataFrame()
)


def determine_pangenome_samples(wildcards):
    return _bpc_manifest.resolve_graph_samples(manifest_list, wildcards.graph, config)


def determine_pangenome_subset(wildcards):
    samples = determine_pangenome_samples(wildcards)
    return metadata.filter(pl.col("Filename").is_in(samples))


# Assembly ingestion & formatting
include: "snakepit/ingest.smk"
include: "snakepit/pangenome.smk"
# Assembly level analyses
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
        "analyses/QC/BPC_metadata_with_QC.All.tsv",
        "analyses/QC/PCA/Bovinae.biSNPs.eigenvec",
        "analyses/satellites/Bovinae.csv",
        "analyses/pangene/Bovinae.ARS_UCD2.0_Ensembl.clustered.tsv",
        "analyses/minigraph/Bovina/L50/mg.gfa",
        "cactus/BosTaurus_n1/graph.gbz",


rule ingest:
    input:
        expand(
            str(ASSEMBLIES_DIR / "{sample}.fa.gz"),
            sample=[r.sample for r in manifest_list],
        ),
        "analyses/ingest/ingest_summary.tsv",


rule qc:
    input:
        "analyses/QC_summary.All.csv",
        "analyses/QC/BPC_metadata_with_QC.All.tsv",
