def _resolve_raw_source(wildcards):
    rec = manifest_records[wildcards.sample]
    if rec.source:
        return rec.source
    for ext in (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna"):
        candidate = RAW_ASSEMBLIES_DIR / f"{wildcards.sample}{ext}"
        if candidate.exists():
            return str(candidate)
    return str(RAW_ASSEMBLIES_DIR / f"{wildcards.sample}.fasta.gz")


def _stage_local_input(wildcards):
    src = _resolve_raw_source(wildcards)
    if re.match(r"^(?:https?|ftp)://|^agc://|^GC[AF]_\d{9}(?:\.\d+)?$", src, re.IGNORECASE):
        return []
    return [src]


_REFERENCE_FASTAS = [
    str(p) for p in config.get("references", {}).values() if str(p).strip()
]

if _REFERENCE_FASTAS:
    rule index_reference:
        input:
            "{ref_fasta}",
        output:
            "{ref_fasta}.fai",
        wildcard_constraints:
            ref_fasta="|".join(re.escape(p) for p in _REFERENCE_FASTAS),
        localrule: True
        shell:
            """
samtools faidx {input}
            """


rule stage_assembly:
    input:
        _stage_local_input,
    output:
        fasta=temp(
            multiext(
                str(STAGED_ASSEMBLIES_DIR / "{sample}.fa.gz"), "", ".fai", ".gzi"
            )
        ),
        ncbi_chrom_map="analyses/ingest/{sample}.ncbi_chrom_map.tsv",
        report="analyses/ingest/{sample}.stage.json",
    params:
        script=workflow.source_path("../scripts/stage_assembly.py"),
        source=_resolve_raw_source,
        checksum=lambda wildcards: manifest_records[wildcards.sample].checksum,
        min_contig_len=config.get("min_contig_len", 0),
        gap_len=config.get("cut_gap_min_len", 10),
    threads: 2
    resources:
        mem_mb_per_cpu=2500,
        runtime="1h",
    shell:
        """
python3 {params.script} \
  --sample {wildcards.sample} \
  --source "{params.source}" \
  --output {output.fasta[0]} \
  --report {output.report} \
  --chrom-map-out {output.ncbi_chrom_map} \
  --checksum "{params.checksum}" \
  --min-contig-len {params.min_contig_len} \
  --gap-len {params.gap_len} \
  --threads {threads}
        """


rule gather_ingest_report:
    input:
        reports=expand(
            "analyses/ingest/{sample}.stage.json",
            sample=[r.sample for r in manifest_list],
        ),
    output:
        tsv="analyses/ingest/ingest_summary.tsv",
    localrule: True
    run:
        cols = [
            "sample",
            "n_contigs",
            "total_bp",
            "n_bases",
            "n_gap_runs",
            "scaffold_n50",
            "scaffold_ng50",
            "contig_n50",
            "contig_ng50",
            "dropped_contigs",
            "dropped_bp",
            "min_len",
            "max_len",
            "gc_frac",
            "source_sha256",
            "source",
        ]
        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output.tsv, "w", encoding="utf-8") as fout:
            fout.write("\t".join(cols) + "\n")
            for fname in input.reports:
                data = json.loads(Path(fname).read_text(encoding="utf-8"))
                fout.write("\t".join(str(data.get(c, "")) for c in cols) + "\n")
