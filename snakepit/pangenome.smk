import polars as pl
from pathlib import PurePath

rule cut_assemblies_at_gaps:
    input:
        fasta = multiext('data/raw_assemblies/{sample}.fasta.gz','','.fai','.gzi')
    output:
        fasta = 'data/contigs/{sample}.fa'
    threads: 1
    resources:
        mem_mb_per_cpu = 1000,
        runtime = '30m'
    shell:
        '''
seqtk cutN -n 0 {input.fasta[0]} |\
sed 's/_RagTag//' > {output}
        '''

# We ideally want a reference with an X, Y, and MT so nothing goes missing...
rule ragtag_scaffold:
    input:
        fasta = rules.cut_assemblies_at_gaps.output['fasta'],
        reference = lambda wildcards: config['references'][wildcards.reference]
    output:
        multiext('analyses/scaffolding/{reference}/{sample}/ragtag.scaffold','.agp','.fasta','.err','.confidence.txt','.stats','.asm.paf','.asm.paf.log')
    params:
        _dir = lambda wildcards, output: PurePath(output[0]).parent,
        mm2_opt = '-x asm20',
        exclude_unplaced = f"^({'|'.join(list(map(str,range(1,30))) + ['X','Y','MT'])})"
    conda: 'RagTag'
    threads: 6
    resources:
        mem_mb_per_cpu = 12000,
        runtime = '2h'
    shell:
        '''
grep -vwP "{params.exclude_unplaced}" {input.reference}.fai > $TMPDIR/unplaced.txt

ragtag.py scaffold {input.reference} {input.fasta} \
  -o {params._dir} \
  --mm2-params "{params.mm2_opt} -t {threads}" \
  -e $TMPDIR/unplaced.txt
        '''

# Uses the layout defined in the BPC spreadsheet
def map_ID_to_filename(sample):
    filename_map = dict(pl.read_csv(config['metadata'])
                       .select(['ID','Filename'])
                       .iter_rows())

    return filename_map[sample]

rule panSN_renaming:
    input:
        metadata = config['metadata'],
        fasta = lambda wildcards: expand(rules.cut_assemblies_at_gaps.output['fasta'],sample=map_ID_to_filename(wildcards.sample)),
        agp = lambda wildcards: expand(rules.ragtag_scaffold.output[0],sample=map_ID_to_filename(wildcards.sample),reference='ARS_UCD2.0')
    output:
        fasta = multiext('data/freeze_1/{sample}.fa.gz','','.fai','.gzi')
    params:
        haplotype = '0' # currently fix as haploid
    threads: 2
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '1h'
    shell:
        '''
seqtk seq -l 0 {input.fasta} |\
awk 'function revcomp(arg) {{o = "";for(i = length(arg); i > 0; i--) {{o = o c[substr(arg, i, 1)]}} return(o)}}; BEGIN {{c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A"}}; {{if (NR==FNR) {{if($5=="W") {{if ($1~/_RagTag/) {{sub("_RagTag","",$1);}}else {{$1="unplaced";}}V[">"$6]=$8;C[$1]++;R[">"$6]=">{wildcards.sample}#{params.haplotype}#"$1"#"C[$1]-1;}}}}else {{if ($1~/^>/) {{printf "%s\\t",R[$1]; Z=V[$1] }} }}else {{if (Z=="+") {{print $1;}} else {{print revcomp($1);}}}}}}}}' {input.agp} - |\
sort -k1,1V |\
tr "\\t" "\\n" |\
bgzip -c -@ 2 - > {output.fasta[0]}

samtools faidx {output.fasta[0]}
        '''

rule agc_create:
    input:
        assemblies = lambda wildcards: expand(rules.panSN_renaming.output['fasta'][0] if wildcards.archive == 'panSN' else 'data/raw_assemblies/{sample}.fasta.gz',sample=pl.read_csv(config['metadata']).get_column('ID').to_list())
    output:
        agc = 'data/freeze_1/{archive}.agc'
    threads: 8
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '4h'
    shell:
        '''
agc create -d -t {threads} {input.assemblies} > {output.agc}
        '''


# something to split fasta output
rule panSN_split:
    input:
        fasta = multiext('data/freeze_1/{sample}.fa.gz','','.fai','.gzi')
    output:
        fasta = expand('data/freeze_1/chromosomes/{{sample}}.{chromosome}.fa.gz{ext}',ext=('','.fai','.gzi'),chromosome=list(map(str,range(1,30))))
    shell:
        '''
for C in {{1..29}}
do
  samtools faidx --continue --write-index -r <(grep "#${{C}}_" {input.fasta[1]} | cut -f 1) -o data/freeze_1/chromosomes/{wildcards.sample}.${{C}}.fa.gz --length 0 {input.fasta[0]}
done
        '''

rule panSN_gather:
    input:
        assemblies = lambda wildcards: expand('data/freeze_1/chromosomes/{sample}.{chromosome}.fa.gz',sample=determine_pangenome_samples(wildcards.graph),allow_missing=True)
    output:
        fasta = multiext('data/freeze_1/{graph}/{chromosome}.fa.gz','','.fai','.gzi')
    shell:
        '''
cat {input.assemblies} > {output.fasta[0]}

samtools faidx {output.fasta[0]}
        '''
