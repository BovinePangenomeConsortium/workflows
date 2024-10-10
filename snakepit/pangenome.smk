# Uses the layout defined in the BPC spreadsheet
def map_ID_to_filename(**wildcards):
    filename_map = dict(pl.read_csv('metdata.csv')
                       .select(['ID','Filename'])
                       .iter_rows())
    return filename_map[wildcards.sample]

rule panSN_renaming:
    input:
        fasta = rules.cut_assemblies_at_gaps.output['fasta'],
        agp = lambda wildcards: expand(rules.ragtag_scaffold.output[0],sample=map_ID_to_filename(wildcards.sample))
    output:
        fasta = multiext('data/freeze_1/{sample}.fa.gz','','.fai','.gzi')
    params:
        haplotype = '0' # currently fix as haploid
    threads: 2
    resources:
        mem_mb_per_cpu = 1500,
        runtime = '1h'
    shell:
        '''
        pigz -dc {input.fasta} |\
        awk 'NR==FNR {{ if($5=="W") {{ sub("_RagTag","",$1); C[$1]++; R[">"$6]=">{wildcards.sample}#{params.haplotype}#"$1"_"C[$1]-1}}; next}} ($1 in R && $1~/^>/) {{$1=R[$1]}}1' {input.agp} - |\
        bgzip -c -@ 2 - > {output.fasta[0]}
        samtools faidx {output.fasta[0]}
        '''

# this implicitly will use the first assembly as the reference to compress against
rule agc_create:
    input:
        assem
        blies = expand(rules.panSN_renaming.output['fasta'],samples=samples)
    output:
        agc = 'data/freeze_1/BPC_freeze_1.agc'
    threads: 8
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '4h'
    shell:
        '''
        agc create -d -t {threads} {input.assemblies} > {output.agc}
        '''

rule mash_triangle:
    input:
        rules.panSN_spec.output
    output:
        'pangenome/tree/{chromosome}.matrix'
    threads: 4
    resources:
        mem_mb = 2500,
        walltime = '1h'
    shell:
        '''
        mash triangle -i -p {threads} {input[0]} > {output}
        '''
