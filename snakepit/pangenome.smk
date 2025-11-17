rule cut_assemblies_at_gaps:
    input:
        fasta=multiext("data/raw_assemblies/{sample}.fasta.gz", "", ".fai", ".gzi"),
    output:
        fasta="data/contigs/{sample}.fa",
    threads: 1
    resources:
        mem_mb_per_cpu=1000,
        runtime="30m",
    shell:
        """
seqtk cutN -n 0 {input.fasta[0]} |\
sed 's/_RagTag//' > {output}
        """


# We ideally want a reference with an X, Y, and MT so nothing goes missing...
rule ragtag_scaffold:
    input:
        fasta=lambda wildcards: (
            "data/raw_assemblies/{sample}.fasta.gz"
            if wildcards.sample in ANNOTATED_GENOMES
            else rules.cut_assemblies_at_gaps.output["fasta"]
        ),
        reference=lambda wildcards: config["references"][wildcards.reference],
    output:
        multiext(
            "analyses/scaffolding/{reference}/{sample}/ragtag.scaffold",
            ".agp",
            ".fasta",
            ".err",
            ".confidence.txt",
            ".stats",
            ".asm.paf",
            ".asm.paf.log",
        ),
    params:
        _dir=lambda wildcards, output: Path(output[0]).parent,
        mm2_opt="--cs -c -x asm10",
        exclude_unplaced=f"^({'|'.join(list(map(str, range(1,30)))+['X', 'Y', 'MT'])})",
        annotated_exclude=lambda wildcards: (
            "-j $TMPDIR/query_exclude.txt"
            if wildcards.sample in ANNOTATED_GENOMES
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
grep -vwP "{params.exclude_unplaced}" {input.reference}.fai | cut -f 1 > $TMPDIR/unplaced.txt
grep -vwP "{params.exclude_unplaced}" {input.fasta}.fai | cut -f 1 > $TMPDIR/query_exclude.txt

ragtag.py scaffold {input.reference} {input.fasta} \
  -o {params._dir} \
  --mm2-params "{params.mm2_opt} -t {threads}" \
  -e $TMPDIR/unplaced.txt {params.annotated_exclude}
        """


# Uses the layout defined in the BPC spreadsheet
def map_ID_to_filename(sample):
    filename_map = dict(
        pl.read_csv(config["metadata"]).select(["ID", "Filename"]).iter_rows()
    )

    return filename_map[sample]


# TODO: need to convert the P lines later in pggb to vg format [start-end] rather than :start-end
def panSN_naming_schema(sample, schema="cactus"):
    haplotype = (
        metadata.filter(pl.col("Filename") == sample)
        .get_column("Haplotype")
        .to_list()[0]
    )
    animal_ID = (
        metadata.filter(pl.col("Filename") == sample)
        .get_column("Animal ID")
        .to_list()[0]
    )
    naming = f'">{animal_ID}#{haplotype}#"$1'
    match schema:
        case "cactus":
            naming += '":"$2"-"$3'
        case "vg":
            naming += '"["$2"-"$3"]"'
        case "_":
            naming += '"_"$2"-"$3'
    return naming


rule panSN_renaming:
    input:
        fasta=rules.ragtag_scaffold.input["fasta"],
        agp=expand(
            rules.ragtag_scaffold.output[0], reference="ARS_UCD2.0", allow_missing=True
        ),
    output:
        fasta=multiext("data/currated_assemblies/{sample}.fa.gz", "", ".fai", ".gzi"),
    params:
        naming_schema=lambda wildcards: panSN_naming_schema(wildcards.sample),
    threads: 2
    resources:
        mem_mb_per_cpu=2500,
        runtime="1h",
    shell:
        """
seqtk seq -l 0 -U {input.fasta} |\
awk 'function revcomp(arg) {{o = "";for(i = length(arg); i > 0; i--) {{o = o c[substr(arg, i, 1)]}} return(o)}}; BEGIN {{c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A"}}; {{if (NR==FNR) {{if($5=="W") {{if ($1~/^[0-9XYM]/&&$1~/_RagTag/) {{sub("_RagTag","",$1);}}else {{$1="unplaced_"NR;}}V[">"$6]=$9;C[$1]++;R[">"$6]={params.naming_schema};}}}}else {{if ($1~/^>/) {{printf "%s\\t",R[$1]; Z=V[$1] }} else {{if (Z=="+") {{print $1;}} else {{print revcomp($1);}}}}}}}}' {input.agp} - |\
sort -k1,1V |\
tr "\\t" "\\n" |\
bgzip -c -@ 2 - > {output.fasta[0]}

samtools faidx {output.fasta[0]}
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
agc create -d -t {threads} {input.assemblies} > {output.agc}
        """


rule panSN_split:
    input:
        fasta=multiext("data/currated_assemblies/{sample}.fa.gz", "", ".fai", ".gzi"),
    output:
        fasta=expand(
            "data/currated_assemblies/chromosomes/{{sample}}.{chromosome}.fa.gz{ext}",
            ext=("", ".fai", ".gzi"),
            chromosome=list(map(str, range(1, 30))),
        ),
    params:
        regex=lambda wildcards: "" if wildcards.sample in ANNOTATED_GENOMES else r"#\d",
        _out=lambda wildcards, output: Path(output["fasta"][0])
        .with_suffix("")
        .with_suffix("")
        .with_suffix(""),
    threads: 1
    resources:
        mem_mb_per_cpu=2500,
        runtime="30m",
    shell:
        """
for C in {{1..29}}
do
  samtools faidx  --write-index -r <(grep -P "#$C:" {input.fasta[1]} | cut -f 1) -o {params._out}.$C.fa.gz --length 0 {input.fasta[0]}
done
        """


# TODO: this still doesn't like getting non-existant chromosomes
rule panSN_split2:
    input:
        fasta=multiext("data/currated_assemblies/{sample}.fa.gz", "", ".fai", ".gzi"),
    output:
        fasta=expand(
            "data/currated_assemblies/chromosomes/{{sample}}.{chromosome}.fa.gz{ext}",
            ext=("", ".fai", ".gzi"),
            chromosome=("X", "Y", "MT"),
        ),
    params:
        regex=lambda wildcards: "" if wildcards.sample in ANNOTATED_GENOMES else r"#\d",
        _out=lambda wildcards, output: Path(output["fasta"][0])
        .with_suffix("")
        .with_suffix("")
        .with_suffix(""),
    threads: 1
    resources:
        mem_mb_per_cpu=2500,
        runtime="30m",
    shell:
        """
for C in X Y MT
do
  #ugly hack to avoid turning off pipefail. Grep returns exit code of 1 if no matches, which we will allow
  {{ grep -P "#$C{params.regex}\\s" {input.fasta[1]} || test $? = 1; }} | cut -f 1 > $TMPDIR/regions.list
  if [ -s $TMPDIR/regions.list ]
  then
    samtools faidx --continue --write-index --region-file $TMPDIR/regions.list -o {params._out}.$C.fa.gz --length 0 {input.fasta[0]}
  else
    touch {params._out}.$C.fa.gz {params._out}.$C.fa.gz.fai {params._out}.$C.fa.gz.gzi
  fi
done
        """


rule panSN_gather:
    input:
        assemblies=expand(
            "data/currated_assemblies/chromosomes/{sample}.{chromosome}.fa.gz",
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        fasta=multiext(
            "data/currated_assemblies/{graph}/{chromosome}.fa.gz", "", ".fai", ".gzi"
        ),
    shell:
        """
cat {input.assemblies} > {output.fasta[0]}

samtools faidx {output.fasta[0]}
        """
