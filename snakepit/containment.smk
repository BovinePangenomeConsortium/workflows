def get_sample_fastq_from_csv(wildcards):
    df = pl.read_csv(config["SRA"])
    return df.rows_by_key("accession",named=True,unique=True)[wildcards.sample]["fastq"].split()

rule seqtk_prepare:
    input:
        fastq=get_sample_fastq_from_csv
    output:
        fasta=('analyses/skope/samples/{sample}.fa')
    threads: 1
    resources:
        mem_mb_per_cpu=5000,
        runtime="4h"
    shell:
        '''
seqtk mergepe {input.fastq} |\
  seqtk seq -A -L 150 -N > {output.fasta}
        '''

rule skope_index:
    input:
        fasta=expand(rules.panSN_renaming.output["fasta"][0],sample=determine_pangenome_samples,allow_missing=True),
    output:
        index="analyses/skope/indexes/{graph}.k{kmer}.s{smer}.idx"
    threads: 4
    resources:
        mem_mb_per_cpu=15000,
        runtime="24h"
    shell:
        '''
ln -s $(realpath {input.fasta}) $TMPDIR
skope index build -t {threads} -o {output.index} -k {wildcards.kmer} -s {wildcards.smer} $TMPDIR
        '''

rule skope_contain:
    input:
        fasta=rules.seqtk_prepare.output["fasta"],
        index=rules.skope_index.output["index"]
    output:
        hits="analyses/skope/containment/{graph}.k{kmer}.s{smer}.{sample}.{mode,group-counts|per-seq}.tsv"
    threads: 4
    resources:
        mem_mb_per_cpu=12500
    shell:
        '''
skope classify -t {threads} --output-mode {wildcards.mode} --names {wildcards.sample} -o {output.hits} {input.index} {input.fasta}
        '''
