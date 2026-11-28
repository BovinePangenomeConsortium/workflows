# Bovine Pangenome Consortium workflows

This is the home for pipelines related to the [Bovine Pangenome Consortium](https://bovinepangenomeconsortium.github.io/).
This is a work in progress, but reach out if you have any questions!

### Using the pipeline

Currently, the pipelines are optimised to run on ETH's cluster Euler.
Tools are assumed to be installed, but eventually conda environments will be supported.

We can test the pipeline itself works (without checking the individual jobs work) with

```
git clone git@github.com:BovinePangenomeConsortium/workflows.git
cd .test
snakemake -s ../Snakefile --configfile config.yaml -n
```

### Adding and formatting new assemblies

Assembly ingestion is driven by the manifest configured in `config.yaml` (e.g., `assemblies.tsv`, or legacy `BPC_metadata.csv`) and declarative graph definitions under `graphs:`.

1. **Register a new assembly** (single sample, paired diploid `--hap1`/`--hap2`, or `--batch incoming.tsv`) using `scripts/add_assembly.py`:
   ```bash
   # Single assembly with explicit animal ID
   python3 scripts/add_assembly.py \
     --manifest assemblies.tsv \
     --source /path/to/asm.fasta.gz \
     --animal-id BOV0123 \
     --haplotype 1 \
     --sex M \
     --species "Bos taurus" \
     --breed Holstein

   # Paired diploid assembly with automated Excel-safe hex ID + private crosswalk
   python3 scripts/add_assembly.py \
     --manifest assemblies.tsv \
     --hap1 /path/to/hap1.fa.gz \
     --hap2 /path/to/hap2.fa.gz \
     --external-id FARM_ANIMAL_42 \
     --id-map private_id_map.tsv \
     --sex M \
     --species "Bos taurus"
   ```
   Omit `--animal-id` (or pass `--animal-id auto`) to mint a collision-free 8-character hex ID (guaranteed to contain at least one `A/B/C/D/F` so spreadsheets/R do not coerce it to a number or scientific notation). `source` accepts local FASTAs (`.fa`/`.fasta`/`.fna`, plain, gzip/bgzip, or zstd), `https://`/`ftp://` URLs, NCBI accessions (`GCA_*`/`GCF_*` via `datasets`), or `agc://archive.agc#sample`.

2. **Choose a `format_mode`** (defaults to `intact`, or `ragtag_preserve` when `reference_annotation=Y`):
   - `intact` (recommended default for Minigraph-Cactus & Pangene): keep all scaffolds intact without gap-cutting or RagTag, map canonical chromosome names (via optional `chrom_map` or `chr*`/`M` normalization), and emit `AnimalID#Hap#Contig` PanSN headers.
   - `already_placed`: skip RagTag and map existing chromosome names with PanSN subrange coordinates (`:start-end` or `[0start-end]`).
   - `ragtag_preserve`: keep existing scaffolds intact and scaffold/orient placed chromosomes with RagTag.
   - `ragtag_cut`: split at gaps (`cut_gap_min_len`, default `10`), scaffold with RagTag, then format to PanSN (useful for PGGB chromosome partitioning).

   Set `min_contig_len` in `config.yaml` (e.g. `10000` or `50000`) to filter tiny unplaced fragments during `stage_assembly` (canonical `MT`/`chrM` contigs are preserved automatically).

3. **Run ingestion and QC** before triggering full pangenome builds:
   ```bash
   snakemake -s Snakefile --configfile config.yaml ingest
   snakemake -s Snakefile --configfile config.yaml qc
   ```
   `qc` writes `analyses/QC_summary.All.csv` and joins it cleanly with the input manifest into `analyses/QC/BPC_metadata_with_QC.All.tsv` with `inferred_sex` and `auto_qc_flags` (leaving your input manifest untouched).

4. **Graph subsets** are defined declaratively in `config.yaml` under `graphs:` (filtering by `exclude_qc_fail`, `species`, `subspecies`, `breed`, `sex`, `tags`, `one_per_animal`, and `include`/`exclude` overrides).

### Details

More details on the motivations and details of the pipeline can be found in the [BPC/docs](https://github.com/BovinePangenomeConsortium/docs) repository.
