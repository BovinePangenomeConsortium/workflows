rule KMC_count:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
    output:
        kmers=multiext(
            "analyses/satellites/{sample}.kmc", ".kmc_pre", ".kmc_suf", ".txt"
        ),
    params:
        prefix=lambda wildcards, output: Path(output[0]).with_suffix(""),
    threads: 2
    resources:
        mem_mb_per_cpu=10000,
        runtime="1h",
    shell:
        """
kmc -k151 -t{threads} -ci20 -cs100000 -fm {input.fasta[0]} {params.prefix} $TMPDIR

kmc_dump {params.prefix} {output.kmers[2]}
        """


# Note: we set `-l 40` since the smallest expected bovine SAT is 46 bp, and to avoid picking up the 23 bp sub-motifs in other SATs
rule SRF:
    input:
        kmers=rules.KMC_count.output["kmers"][2],
    output:
        satellites="analyses/satellites/{sample}.srf.fa",
    resources:
        mem_mb_per_cpu=2500,
        runtime="30m",
    shell:
        """
        srf -p {wildcards.sample} -l 40 -c 100 {input.kmers} > {output.satellites}
        """


rule minimap2_srf:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
        satellites=rules.SRF.output["satellites"],
    output:
        paf="analyses/satellites/{sample}.srf-aln.paf",
    threads: 6
    resources:
        mem_mb_per_cpu=2500,
        runtime="4h",
    shell:
        """
minimap2 -t {threads} -c -N1000000 -f1000 -r100,100 <(srfutils.js enlong {input.satellites}) {input.fasta[0]} > {output.paf}
        """


rule srfutils:
    input:
        paf=rules.minimap2_srf.output["paf"],
        fai=rules.panSN_renaming.output["fasta"][1],
    output:
        bed="analyses/satellites/{sample}.srf-aln.bed",
        abun="analyses/satellites/{sample}.srf-aln.abun",
    threads: 1
    resources:
        mem_mb_per_cpu=5000,
        runtime="1h",
    shell:
        """
        srfutils.js paf2bed {input.paf} > {output.bed}

        gsize=$(awk '{{s+=$2}} END {{print s}}' {input.fai})
        srfutils.js bed2abun -g "$gsize" {output.bed} > {output.abun}
        """


rule srf_gather:
    input:
        abun=expand(
            rules.srfutils.output["abun"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
        satellites=expand(
            rules.SRF.output["satellites"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
        fai=expand(
            rules.panSN_renaming.output["fasta"][1],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        csv="analyses/satellites/{graph}.csv",
    localrule: True
    run:
        srf_pattern = re.compile(
            r">(?P<sample>[^\s#]+)#circ(?P<contig>\d+)-(?P<length>\d+) "
            r"min=(?P<min>\d+),max=(?P<max>\d+),avg=(?P<avg>\d+)"
        )
        sat_rows = []
        for fname in input.satellites:
            with open(fname) as fin:
                for line in fin:
                    if line.startswith(">"):
                        d = srf_pattern.match(line.rstrip()).groupdict()
                        sat_rows.append(
                            {
                                "sample": d["sample"],
                                "contig": d["contig"],
                                "length": int(d["length"]),
                                "min": int(d["min"]),
                                "max": int(d["max"]),
                                "avg": int(d["avg"]),
                            }
                        )
                    elif sat_rows:
                        sat_rows[-1]["fasta"] = line.rstrip()
        satellites = pl.DataFrame(
            sat_rows,
            schema={
                "sample": pl.String,
                "contig": pl.String,
                "length": pl.Int64,
                "min": pl.Int64,
                "max": pl.Int64,
                "avg": pl.Int64,
                "fasta": pl.String,
            },
        )

        abun_pattern = re.compile(
            r"(?P<sample>[^\s#]+)#circ(?P<contig>\d+)-(?P<length>\d+)\t"
            r"(?P<bp>\d+)\t(?P<mean_divergence>\S+)\t(?P<filtered>\S+)\t(?P<unfiltered>\S+)"
        )
        abun_rows = []
        for fname in input.abun:
            with open(fname) as fin:
                for line in fin:
                    m = abun_pattern.match(line.rstrip())
                    if m:
                        d = m.groupdict()
                        div_str = d["mean_divergence"]
                        abun_rows.append(
                            {
                                "sample": d["sample"],
                                "contig": d["contig"],
                                "length": int(d["length"]),
                                "bp": int(d["bp"]),
                                "mean_divergence": (
                                    -1.0 if div_str == "NaN" else float(div_str)
                                ),
                                "filtered": float(d["filtered"]),
                                "unfiltered": float(d["unfiltered"]),
                            }
                        )
        abundances = pl.DataFrame(
            abun_rows,
            schema={
                "sample": pl.String,
                "contig": pl.String,
                "length": pl.Int64,
                "bp": pl.Int64,
                "mean_divergence": pl.Float64,
                "filtered": pl.Float64,
                "unfiltered": pl.Float64,
            },
        )

        fai_rows = []
        for fname in input.fai:
            sample = Path(fname).name.removesuffix(".fa.gz.fai")
            total_bp = 0
            with open(fname) as fin:
                for line in fin:
                    if line.strip():
                        total_bp += int(line.split("\t")[1])
            fai_rows.append({"sample": sample, "assembly_bp": total_bp})
        fai_df = pl.DataFrame(
            fai_rows,
            schema={"sample": pl.String, "assembly_bp": pl.Int64},
        )

        (
            satellites.join(abundances, on=["sample", "contig", "length"])
            .join(fai_df, on="sample", how="left")
            .with_columns(
                (pl.col("bp") / pl.col("assembly_bp")).alias("genome_frac"),
                (pl.col("bp") / pl.col("bp").sum().over("sample")).alias("sat_frac"),
            )
            .write_csv(output["csv"])
        )


rule srf_build_catalog:
    input:
        csv=rules.srf_gather.output["csv"],
        refs=workflow.source_path("../resources/bovine_satellite_refs.fa"),
    output:
        catalog="analyses/satellites/{graph}.catalog.fa",
        membership="analyses/satellites/{graph}.membership.tsv",
        summary="analyses/satellites/{graph}.catalog_summary.tsv",
    params:
        script=workflow.source_path("../scripts/cluster_srf_motifs.py"),
        prefix="analyses/satellites/{graph}",
    localrule: True
    shell:
        """
python3 {params.script} \
  --csv {input.csv} \
  --refs {input.refs} \
  --out-prefix {params.prefix}
        """


rule srf_chrom_summary:
    input:
        bed=rules.srfutils.output["bed"],
        fai=rules.panSN_renaming.output["fasta"][1],
        membership=rules.srf_build_catalog.output["membership"],
    output:
        tsv="analyses/satellites/chrom_summary/{graph}/{sample}.tsv",
    params:
        script=workflow.source_path("../scripts/summarize_srf_bed.py"),
    threads: 1
    resources:
        mem_mb_per_cpu=2500,
        runtime="30m",
    shell:
        """
python3 {params.script} \
  --bed {input.bed} \
  --fai {input.fai} \
  --sample {wildcards.sample} \
  --membership {input.membership} \
  --output {output.tsv}
        """


rule srf_chrom_summary_gather:
    input:
        tsv=expand(
            rules.srf_chrom_summary.output["tsv"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        tsv="analyses/satellites/{graph}.chrom_summary.tsv",
    localrule: True
    shell:
        """
awk 'FNR>1||NR==1' {input.tsv} > {output.tsv}
        """


rule minimap2_srf_catalog:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
        catalog=rules.srf_build_catalog.output["catalog"],
    output:
        paf="analyses/satellites/unified/{graph}/{sample}.paf",
    threads: 6
    resources:
        mem_mb_per_cpu=2500,
        runtime="4h",
    shell:
        """
minimap2 -t {threads} -c -N1000000 -f1000 -r100,100 <(srfutils.js enlong {input.catalog}) {input.fasta[0]} > {output.paf}
        """


rule srfutils_catalog:
    input:
        paf=rules.minimap2_srf_catalog.output["paf"],
        fai=rules.panSN_renaming.output["fasta"][1],
    output:
        bed="analyses/satellites/unified/{graph}/{sample}.bed",
        abun="analyses/satellites/unified/{graph}/{sample}.abun",
    threads: 1
    resources:
        mem_mb_per_cpu=5000,
        runtime="1h",
    shell:
        """
        srfutils.js paf2bed {input.paf} > {output.bed}

        gsize=$(awk '{{s+=$2}} END {{print s}}' {input.fai})
        srfutils.js bed2abun -g "$gsize" {output.bed} > {output.abun}
        """


rule srf_catalog_chrom_summary:
    input:
        bed=rules.srfutils_catalog.output["bed"],
        fai=rules.panSN_renaming.output["fasta"][1],
        catalog=rules.srf_build_catalog.output["catalog"],
    output:
        tsv="analyses/satellites/unified/{graph}/{sample}.chrom_summary.tsv",
    params:
        script=workflow.source_path("../scripts/summarize_srf_bed.py"),
    threads: 1
    resources:
        mem_mb_per_cpu=2500,
        runtime="30m",
    shell:
        """
python3 {params.script} \
  --bed {input.bed} \
  --fai {input.fai} \
  --sample {wildcards.sample} \
  --catalog {input.catalog} \
  --output {output.tsv}
        """


rule srf_catalog_gather:
    input:
        abun=expand(
            rules.srfutils_catalog.output["abun"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
        chrom_sum=expand(
            rules.srf_catalog_chrom_summary.output["tsv"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
        summary=rules.srf_build_catalog.output["summary"],
        fai=expand(
            rules.panSN_renaming.output["fasta"][1],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        abun_tsv="analyses/satellites/{graph}.unified_abun.tsv",
        chrom_tsv="analyses/satellites/{graph}.unified_chrom_summary.tsv",
    localrule: True
    run:
        abun_rows = []
        for fname in input.abun:
            sample = Path(fname).stem
            with open(fname) as fin:
                for line in fin:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 5:
                        continue
                    abun_rows.append(
                        {
                            "sample": sample,
                            "catalog_id": fields[0],
                            "bp": int(fields[1]),
                            "mean_divergence": (
                                -1.0 if fields[2] == "NaN" else float(fields[2])
                            ),
                            "filtered": float(fields[3]),
                            "unfiltered": float(fields[4]),
                        }
                    )
        abun_df = pl.DataFrame(
            abun_rows,
            schema={
                "sample": pl.String,
                "catalog_id": pl.String,
                "bp": pl.Int64,
                "mean_divergence": pl.Float64,
                "filtered": pl.Float64,
                "unfiltered": pl.Float64,
            },
        )

        fai_rows = []
        for fname in input.fai:
            sample = Path(fname).name.removesuffix(".fa.gz.fai")
            total_bp = 0
            with open(fname) as fin:
                for line in fin:
                    if line.strip():
                        total_bp += int(line.split("\t")[1])
            fai_rows.append({"sample": sample, "assembly_bp": total_bp})
        fai_df = pl.DataFrame(
            fai_rows,
            schema={"sample": pl.String, "assembly_bp": pl.Int64},
        )

        cat_summary = pl.read_csv(input.summary, separator="\t").select(
            ["catalog_id", "family_id", "subfamily_id", "ref_name", "ref_acc", "rep_length"]
        )

        (
            abun_df.join(cat_summary, on="catalog_id", how="left")
            .join(fai_df, on="sample", how="left")
            .with_columns(
                (pl.col("bp") / pl.col("assembly_bp")).alias("genome_frac"),
                (pl.col("bp") / pl.col("bp").sum().over("sample")).alias("sat_frac"),
            )
            .write_csv(output["abun_tsv"], separator="\t")
        )

        with open(output["chrom_tsv"], "w") as fout:
            for idx, fname in enumerate(input.chrom_sum):
                with open(fname) as fin:
                    for line_num, line in enumerate(fin):
                        if line_num == 0 and idx > 0:
                            continue
                        fout.write(line)


rule btsat1_genotype:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
    output:
        alleles="analyses/satellites/btsat1/{sample}.alleles.tsv",
        spacing="analyses/satellites/btsat1/{sample}.spacing.tsv",
    params:
        script=workflow.source_path("../scripts/genotype_btsat1_alleles.py"),
    threads: 1
    resources:
        mem_mb_per_cpu=2500,
        runtime="1h",
    shell:
        """
python3 {params.script} \
  --fasta {input.fasta[0]} \
  --sample {wildcards.sample} \
  --out-alleles {output.alleles} \
  --out-spacing {output.spacing}
        """


rule btsat1_genotype_gather:
    input:
        alleles=expand(
            rules.btsat1_genotype.output["alleles"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
        spacing=expand(
            rules.btsat1_genotype.output["spacing"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        alleles="analyses/satellites/{graph}.btsat1_alleles.tsv",
        spacing="analyses/satellites/{graph}.btsat1_spacing.tsv",
    localrule: True
    shell:
        """
awk 'FNR>1||NR==1' {input.alleles} > {output.alleles}
awk 'FNR>1||NR==1' {input.spacing} > {output.spacing}
        """


rule kmc_check_btsat1:
    input:
        kmers=rules.KMC_count.output["kmers"][2],
    output:
        tsv="analyses/satellites/btsat1/{sample}.kmc_check.tsv",
    params:
        script=workflow.source_path("../scripts/check_kmc_btsat1.py"),
    threads: 1
    resources:
        mem_mb_per_cpu=2500,
        runtime="30m",
    shell:
        """
python3 {params.script} \
  --kmc {input.kmers} \
  --sample {wildcards.sample} \
  --output {output.tsv}
        """


rule kmc_check_btsat1_gather:
    input:
        tsv=expand(
            rules.kmc_check_btsat1.output["tsv"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        tsv="analyses/satellites/{graph}.btsat1_kmc_check.tsv",
    localrule: True
    shell:
        """
awk 'FNR>1||NR==1' {input.tsv} > {output.tsv}
        """
