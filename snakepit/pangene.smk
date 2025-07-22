rule miniprot_index:
    input:
        fasta = rules.panSN_renaming.output['fasta']
    output:
        index = 'analyses/pangene/{sample}.mpi'
    threads: 4
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '30m'
    shell:
        '''
miniprot -t {threads} -d {output.index} {input.fasta}
        '''

rule download_annotated_peptides:
    output:
        fasta = 'analyses/pangene/peptides/{peptides}.pep.fa.gz'
    params:
        url = lambda wildcards: config['peptides'][wildcards.peptides]
    localrule: True
    shell:
        '''
wget -O {output.fasta} {params.url} 
        '''

rule minisplice_model:
    output:
        multiext('vi2-7k.kan','','.cali','.log')
    localrule: True
    shell:
        '''
wget -O- https://zenodo.org/records/15931054/files/vi2-7k.tgz | tar zxf - 
        '''

rule minisplice_predict:
    input:
        fasta = rules.panSN_renaming.output['fasta'],
        cali = rules.minisplice_model.output
    output:
        splicing = 'analyses/pangene/{sample}.tsv'
    threads: 8
    resources:
        mem_mb_per_cpu = 1000
    shell:
        '''
minisplice predict -t {threads} -c {input.cali[1]} {input.cali[0]} {input.fasta[0]} > {output.splicing}
        '''
        
rule miniprot_align:
    input:
        fasta = rules.miniprot_index.output['index'] if config.get('miniprot_index',False) else rules.panSN_renaming.output['fasta'],
        splicing = rules.minisplice_predict.output['splicing'],
        peptides = rules.download_annotated_peptides.output['fasta'] 
    output:
        paf = 'analyses/pangene/{sample}.{peptides}.paf.gz'
    threads: 8
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '2h'
    shell:
        '''
miniprot -t {threads} --outs=0.97 -I -u --spsc={input.splicing} {input.fasta} {input.peptides} |\
pigz -p {threads} -c > {output.paf}
        '''

rule pangene:
    input:
        paf = expand(rules.miniprot_align.output['paf'],sample=determine_pangenome_samples,allow_missing=True)
    output:
        gfa = 'analyses/pangene/{graph}.{peptides}.gfa'
    threads: 1
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '1h'
    shell:
        '''
for F in {input.paf}
do
  zgrep -vhP "#(X|Y|MT|unplaced)#" $F | pigz -c > $TMPDIR/$(basename $F)
done
pangene $TMPDIR/*paf.gz > {output.gfa}
        '''

rule pangene_matrix:
    input:
        gfa = rules.pangene.output['gfa']
    output:
        tsv = 'analyses/pangene/{graph}.{peptides}.tsv'
    localrule: True
    shell:
        '''
pangene.js gfa2matrix -c {input.gfa} > {output.tsv}
        '''

rule pangene_call:
    input:
        gfa = rules.pangene.output['gfa']
    output:
        'analyses/pangene/{graph}.{peptides}.call'
    localrule: True
    shell:
        '''
pangene.js call -p -s {input.gfa} > {output}
        '''
