# We want to rename RefSeq but not Ensembl to match 1-30, X, Y, MT naming
rule download_annotated_peptides:
    output:
        multiext(
            "analyses/pangene/peptides/{peptides}",
            peptides=".pep.faa.gz",
            GTF=".gtf.gz",
        ),
    params:
        url_peptides=lambda wildcards: config["peptides"][wildcards.peptides][
            "proteins"
        ],
        url_GTF=lambda wildcards: config["peptides"][wildcards.peptides]["GTF"],
        url_chromalias="https://hgdownload.soe.ucsc.edu/hubs/GCF/002/263/795/GCF_002263795.3/GCF_002263795.3.chromAlias.txt",
    localrule: True
    shell:
        """
wget --quiet -O {wildcards.peptides}.TMP.fa.gz {params.url_peptides}
wget --quiet -O {output.GTF} {params.url_GTF}

if [[ "{wildcards.peptides}" == *"RefSeq"* ]]; then
    wget {params.url_chromalias}
    zcat {wildcards.peptides}.TMP.fa.gz | awk -F'\\t' 'NR==FNR{{if($0!~/^#/){{key=$1;val=$4;if(val=="")val=key;a[key]=val}}next}} !/^#/{{if($1 in a)$1=a[$1]}}1' GCF_002263795.3.chromAlias.txt - | bgzip -c > {output.peptides}
    rm GCF_002263795.3.chromAlias.txt {wildcards.peptides}.TMP.fa.gz
elif [[ "{wildcards.peptides}" == *"Ensembl"* ]]; then
  mv {wildcards.peptides}.TMP.fa.gz {output.peptides}
fi
        """


# No observed difference in protein_coding transcripts for Ensembl, and no record for RefSeq, so don't use `-e`
# Need to use pangene.js fork
rule pangene_getaa:
    input:
        peptides=rules.download_annotated_peptides.output["peptides"],
        GTF=rules.download_annotated_peptides.output["GTF"],
    output:
        multiext(
            "analyses/pangene/peptides/{peptides}.pangene",
            peptides=".faa.gz",
            log=".log",
        ),
    params:
        canonical=lambda wildcards: (
            "-c" if wildcards.peptides.endswith("Ensembl") else ""
        ),
    localrule: True
    shell:
        """
pangene.js getaa {params.canonical} {input.GTF} {input.peptides} 2> {output.log} | bgzip -c > {output.peptides}
        """


rule cd_hit_collapse_genes:
    input:
        peptides=rules.pangene_getaa.output["peptides"],
    output:
        multiext(
            "analyses/pangene/peptides/{peptides}.pep.clustered.faa",
            peptides="",
            clusters=".clstr",
        ),
    threads: 2
    resources:
        mem_mb_per_cpu=2500,
        runtime="30m",
    shell:
        """
cd-hit -i {input.peptides} -o {output.peptides} -c 0.9 -aS 0.8 -n 5 -T {threads} -g
        """


rule minisplice_model:
    output:
        multiext("analyses/pangene/vi2-7k.kan", "", ".cali", ".log"),
    params:
        _dir=lambda wildcards, output: Path(output[0]).parent,
        url="https://zenodo.org/records/15931054/files/vi2-7k.tgz",
    localrule: True
    shell:
        """
mkdir -p {params._dir}

wget -O- {params.url} | tar -xz -C {params._dir}
        """


rule minisplice_predict:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
        cali=rules.minisplice_model.output,
    output:
        splicing="analyses/pangene/{sample}.minisplice.tsv",
    threads: 8
    resources:
        mem_mb_per_cpu=1000,
    shell:
        """
minisplice predict -t {threads} -c {input.cali[1]} {input.cali[0]} {input.fasta[0]} > {output.splicing}
        """


rule miniprot_index:
    input:
        fasta=rules.panSN_renaming.output["fasta"],
    output:
        index="analyses/pangene/{sample}.mpi",
    threads: 4
    resources:
        mem_mb_per_cpu=10000,
        runtime="30m",
    shell:
        """
miniprot -t {threads} -d {output.index} {input.fasta}
        """


rule miniprot_align:
    input:
        fasta=(
            rules.miniprot_index.output["index"]
            if config.get("miniprot_index", False)
            else rules.panSN_renaming.output["fasta"]
        ),
        splicing=rules.minisplice_predict.output["splicing"],
        peptides=lambda wildcards: (
            rules.cd_hit_collapse_genes.output["peptides"]
            if wildcards.clustered == "clustered"
            else rules.pangene_getaa.output["peptides"]
        ),
    output:
        paf="analyses/pangene/{sample}.{peptides}.{clustered}.paf.gz",
    threads: 8
    resources:
        mem_mb_per_cpu=5000,
        runtime="2h",
    shell:
        """
miniprot -t {threads} --outs=0.97 -I -u --spsc={input.splicing} {input.fasta} {input.peptides} |\
pigz -p {threads} -c > {output.paf}
        """


# note: uses custom ASLeonard/pangene fork with panSN naming scheme support
rule pangene:
    input:
        paf=expand(
            rules.miniprot_align.output["paf"],
            sample=determine_pangenome_samples,
            allow_missing=True,
        ),
    output:
        gfa="analyses/pangene/{graph}.{peptides}.{clustered}.gfa",
    threads: 1
    resources:
        mem_mb_per_cpu=5000,
        runtime="1h",
    shell:
        """
pangene -N {input.paf} > {output.gfa}
        """


rule pangene_matrix:
    input:
        gfa=rules.pangene.output["gfa"],
        clusters=lambda wildcards: (
            rules.cd_hit_collapse_genes.output["clusters"]
            if wildcards.clustered == "clustered"
            else []
        ),
    output:
        tsv="analyses/pangene/{graph}.{peptides}.{clustered}.tsv",
    params:
        collapse_paralogs=lambda wildcards, input: (
            f"-d {input.clusters}" if wildcards.clustered == "clustered" else ""
        ),
    localrule: True
    shell:
        """
pangene.js gfa2matrix -c {params.collapse_paralogs} {input.gfa} > {output.tsv}
        """


rule pangene_call:
    input:
        gfa=rules.pangene.output["gfa"],
    output:
        "analyses/pangene/{graph}.{peptides}.{clustered}.call",
    localrule: True
    shell:
        """
pangene.js call -p -s {input.gfa} > {output}
        """
