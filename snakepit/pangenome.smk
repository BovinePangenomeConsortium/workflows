rule cut_assemblies_at_gaps:
    input:
        fasta=rules.stage_assembly.output["fasta"][0],
    output:
        fasta=temp("data/contigs/{sample}.fa"),
    params:
        min_gap=config.get("cut_gap_min_len", 10),
    threads: 1
    resources:
        mem_mb_per_cpu=1000,
        runtime="30m",
    shell:
        """
seqtk cutN -n {params.min_gap} {input.fasta} > {output.fasta}
        """


def _ragtag_inputs(wildcards):
    rec = manifest_records[wildcards.sample]
    ref_path = config["references"][wildcards.reference]
    ref_fai = f"{ref_path}.fai"
    if rec.format_mode in ("ragtag_preserve", "already_placed", "intact"):
        return {
            "fasta": rules.stage_assembly.output["fasta"][0],
            "fai": rules.stage_assembly.output["fasta"][1],
            "reference": ref_path,
            "reference_fai": ref_fai,
        }
    return {
        "fasta": rules.cut_assemblies_at_gaps.output["fasta"],
        "reference": ref_path,
        "reference_fai": ref_fai,
    }


# We ideally want a reference with an X, Y, and MT so nothing goes missing...
rule ragtag_scaffold:
    input:
        unpack(_ragtag_inputs),
    output:
        agp="analyses/scaffolding/{reference}/{sample}/ragtag.scaffold.agp",
        fasta=temp("analyses/scaffolding/{reference}/{sample}/ragtag.scaffold.fasta"),
        err="analyses/scaffolding/{reference}/{sample}/ragtag.scaffold.err",
        confidence="analyses/scaffolding/{reference}/{sample}/ragtag.scaffold.confidence.txt",
        stats="analyses/scaffolding/{reference}/{sample}/ragtag.scaffold.stats",
        paf="analyses/scaffolding/{reference}/{sample}/ragtag.scaffold.asm.paf",
        log="analyses/scaffolding/{reference}/{sample}/ragtag.scaffold.asm.paf.log",
    params:
        _dir=lambda wildcards, output: Path(output.agp).parent,
        mm2_opt=config.get("ragtag_mm2_opt", "--cs -c -x asm10"),
        exclude_unplaced=f"^({'|'.join(ALL_CHROMOSOME)})\\t",
        annotated_exclude=lambda wildcards: (
            "-j $TMPDIR/query_exclude.txt"
            if manifest_records[wildcards.sample].format_mode
            in ("ragtag_preserve", "already_placed", "intact")
            else ""
        ),
    conda:
        "RagTag"
    threads: 6
    resources:
        mem_mb_per_cpu=10000,
        runtime="4h",
    shell:
        """
{{ grep -vP "{params.exclude_unplaced}" {input.reference_fai} || test $? = 1; }} | cut -f 1 > $TMPDIR/unplaced.txt

if [ -n "{params.annotated_exclude}" ];
then
  {{ grep -vP "{params.exclude_unplaced}" {input.fai} || test $? = 1; }} | cut -f 1 > $TMPDIR/query_exclude.txt
fi

ragtag.py scaffold {input.reference} {input.fasta} \
  -o {params._dir} \
  --mm2-params "{params.mm2_opt} -t {threads}" \
  -e $TMPDIR/unplaced.txt {params.annotated_exclude}
        """


def _pansn_renaming_inputs(wildcards):
    rec = manifest_records[wildcards.sample]
    if rec.format_mode in ("intact", "already_placed"):
        inputs = {
            "fasta": rules.stage_assembly.output["fasta"][0],
            "fai": rules.stage_assembly.output["fasta"][1],
        }
        if rec.chrom_map:
            inputs["chrom_map"] = rec.chrom_map
        else:
            inputs["ncbi_chrom_map"] = rules.stage_assembly.output["ncbi_chrom_map"]
        return inputs
    if rec.format_mode == "ragtag_preserve":
        return {
            "fasta": rules.stage_assembly.output["fasta"][0],
            "fai": rules.stage_assembly.output["fasta"][1],
            "agp": expand(
                rules.ragtag_scaffold.output["agp"],
                reference=REFERENCE_ID,
                allow_missing=True,
            )[0],
        }
    return {
        "fasta": rules.cut_assemblies_at_gaps.output["fasta"],
        "agp": expand(
            rules.ragtag_scaffold.output["agp"],
            reference=REFERENCE_ID,
            allow_missing=True,
        )[0],
    }


def _pansn_mode_args(wildcards, input):
    rec = manifest_records[wildcards.sample]
    if rec.format_mode in ("intact", "already_placed"):
        args = "--intact" if rec.format_mode == "intact" else "--already-placed"
        cmap = getattr(input, "chrom_map", "")
        if cmap:
            return f"{args} --chrom-map {cmap}"
        ncbi_cmap = getattr(input, "ncbi_chrom_map", "")
        if ncbi_cmap:
            return f"{args} --optional-chrom-map {ncbi_cmap}"
        return args
    return f"--agp {input.agp}"


rule panSN_renaming:
    input:
        unpack(_pansn_renaming_inputs),
    output:
        fasta=multiext(
            str(ASSEMBLIES_DIR / "{sample}.fa.gz"),
            "",
            ".fai",
            ".gzi",
        ),
        report="analyses/ingest/{sample}.pansn.json",
    params:
        script=workflow.source_path("../scripts/format_pansn.py"),
        animal_id=lambda wildcards: manifest_records[wildcards.sample].animal_id,
        haplotype=lambda wildcards: manifest_records[wildcards.sample].haplotype,
        schema=config.get("pansn_schema", "cactus"),
        chroms=",".join(ALL_CHROMOSOME),
        mode_args=_pansn_mode_args,
    threads: 2
    resources:
        mem_mb_per_cpu=2500,
        runtime="1h",
    shell:
        """
python3 {params.script} format \
  --fasta {input.fasta} \
  --animal-id {params.animal_id} \
  --haplotype {params.haplotype} \
  --output {output.fasta[0]} \
  {params.mode_args} \
  --schema {params.schema} \
  --chroms "{params.chroms}" \
  --threads {threads} \
  --report {output.report}
        """


rule agc_create:
    input:
        assemblies=expand(
            rules.panSN_renaming.output["fasta"][0], sample=determine_pangenome_samples
        ),
    output:
        agc="data/{graph}.agc",
    threads: 8
    resources:
        mem_mb_per_cpu=10000,
        runtime="4h",
    shell:
        """
agc create -d -t {threads} -o {output.agc} {input.assemblies}
        """


rule panSN_split:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
    output:
        fasta=expand(
            ASSEMBLIES_DIR / "chromosomes/{{sample}}.{chromosome}.fa.gz",
            chromosome=ALL_CHROMOSOME,
        ),
        fai=expand(
            ASSEMBLIES_DIR / "chromosomes/{{sample}}.{chromosome}.fa.gz.fai",
            chromosome=ALL_CHROMOSOME,
        ),
        gzi=expand(
            ASSEMBLIES_DIR / "chromosomes/{{sample}}.{chromosome}.fa.gz.gzi",
            chromosome=ALL_CHROMOSOME,
        ),
    params:
        script=workflow.source_path("../scripts/format_pansn.py"),
        out_dir=str(ASSEMBLIES_DIR / "chromosomes"),
        chroms=",".join(ALL_CHROMOSOME),
    threads: 1
    resources:
        mem_mb_per_cpu=2500,
        runtime="30m",
    shell:
        """
python3 {params.script} split \
  --fasta {input.fasta[0]} \
  --out-dir {params.out_dir} \
  --sample {wildcards.sample} \
  --chroms "{params.chroms}" \
  --threads {threads}
        """


rule panSN_gather:
    input:
        assemblies=expand(
            ASSEMBLIES_DIR / "chromosomes/{sample}.{chromosome}.fa.gz",
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        fasta=multiext(
            str(ASSEMBLIES_DIR / "{graph}/{chromosome}.fa.gz"),
            "",
            ".fai",
            ".gzi",
        ),
    shell:
        """
cat {input.assemblies} > {output.fasta[0]}

samtools faidx {output.fasta[0]}
        """
