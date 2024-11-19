

ruleorder: split_approx_mappings_in_chunks > wfmash

wildcard_constraints:
    mode = r'mapping|alignment',
    chunk = r'\.\d+|'

rule wfmash_index:
    input:
        fasta = multiext('data/freeze_1/{graph}/{chromosome}.fa.gz','','.fai','.gzi')
    output:
        index = 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.index.mm3'
    params:
        block_length = lambda wildcards: int(wildcards.segment_length) * 3
    threads: 4
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '1h'
    shell:
        '''
wfmash \
-s {wildcards.segment_length} \
-l {params.block_length} \
-c 2000 \
-P 50k \
-p {wildcards.p} \
-n 1 \
-k 15 \
-t {threads} \
--tmp-base $TMPDIR \
--write-index {output.index} \
{input.fasta[0]}
        '''

rule wfmash:
    input:
        fasta = multiext('data/freeze_1/{graph}/{chromosome}.fa.gz','','.fai','.gzi'),
        index = rules.wfmash_index.output['index'],
        mapping = lambda wildcards: 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.mapping{chunk}.paf' if wildcards.mode == 'alignment' else []
    output:
        paf = 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.{mode}{chunk}.paf'
    params:
        block_length = lambda wildcards: int(wildcards.segment_length) * 3,
        mapping = lambda wildcards, input: f'--align-paf {input.mapping}' if wildcards.mode == 'alignment' else '--approx-mapping'
    threads: lambda wildcards: 12 if wildcards.mode == 'mapping' else 16
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '4h'
    shell:
        '''
wfmash \
-s {wildcards.segment_length} \
-l {params.block_length} \
-c 2000 \
-P 50k \
-p {wildcards.p} \
-n 1 \
-k 15 \
-t {threads} \
--tmp-base $TMPDIR \
--read-index {input.index} \
{params.mapping} \
{input.fasta[0]} \
--lower-triangular \
> {output.paf}
        '''

#weighted random round robin style reimplementation of https://github.com/waveygang/wfmash/blob/main/scripts/split_approx_mappings_in_chunks.py
rule split_approx_mappings_in_chunks:
    input:
        mapping = expand(rules.wfmash.output['paf'],mode='mapping',chunk='',allow_missing=True)
    output:
        paf = expand('analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.mapping.{chunk}.paf',chunk=range(config.get('wfmash_chunks',1)),allow_missing=True)
    params:
        batch_size = 1000 #re-weight the random choices every N mappings
    threads: 1
    resources:
        mem_mb_per_cpu = 1000,
        runtime = '15m'
    run:
        import numpy as np
        rng = np.random.default_rng()
        weights = np.ones(shape=config.get('wfmash_chunks',1))

        fds = [open(fname, 'w') for fname in output['paf']]

        def calculate_p(weights):
            p = 1-weights/weights.sum()
            return p/p.sum()

        indexes = None
        with open(input['mapping'][0]) as fin:
            for rank, line in enumerate(fin):
                if rank % params.batch_size == 0:
                    indexes = rng.choice(config.get('wfmash_chunks',1),params.batch_size,p=calculate_p(weights))
                idx = indexes[rank % params.batch_size]
                fds[idx].write(line)

                _, _, query_start, query_end, _, _, _, target_start, target_end, _, _, _, estimated_identity = line.strip().split('\t')[:13]

                # increase weights by length * divergence
                weights[idx] += max(int(query_end) - int(query_start), int(target_end) - int(target_start)) * ( 100 - float(estimated_identity.split('id:f:')[1]))

        for f in fds:
            f.close()

rule wfmash_concat:
    input:
        expand(rules.wfmash.output['paf'],mode='alignment.',chunk=range(config.get('wfmash_chunks',1)),allow_missing=True)
    output:
        paf = 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.alignment.concat.paf'
    localrule: True
    shell:
        '''
cat {input} > {output}
        '''

rule seqwish:
    input:
        fasta = multiext('data/freeze_1/{graph}/{chromosome}.fa.gz','','.fai','.gzi'),
        alignment = expand(rules.wfmash_concat.output['paf'],allow_missing=True) if config.get('wfmash_chunks',1) > 1 else expand(rules.wfmash.output['paf'],mode='alignment',chunk='',allow_missing=True)
    output:
        gfa = 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.seqwish.gfa'
    threads: 12
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '4h'
    shell:
        '''
seqwish \
-s {input.fasta[0]} \
-p {input.alignment} \
-k {wildcards.k} \
-f 0 \
-g {output.gfa} \
-B 10000000 \
-t {threads} \
--temp-dir $TMPDIR
        '''

def POA_params(wildcards):
    match wildcards.POA:
        case 'asm5':
            return "1,19,39,3,81,1"
        case 'asm10':
            return "1,9,16,2,41,1"
        case 'asm15':
            return "1,7,11,2,33,1"
        case 'asm20':
            return "1,4,6,2,26,1"

rule smoothxg:
    input:
        fasta = multiext('data/freeze_1/{graph}/{chromosome}.fa.gz','','.fai','.gzi'),
        gfa = rules.seqwish.output['gfa']
    output:
        gfa = 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.smoothxg.og'
    params:
        block_id_min = lambda wildcards: round(float(wildcards.p) / 4,4),
        n_haps = lambda wildcards, input: sum(1 for _ in open(input.fasta[1])),
        POA_pad_depth = lambda wildcards, input: 100 * sum(1 for _ in open(input.fasta[1])),
        POA_lengths = '700,900,1100',
        POA_params = POA_params
    threads: 4
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '24h'
    shell:
        '''
smoothxg \
-t {threads} \
-T {threads} \
-g {input.gfa} \
-r {params.n_haps} \
--base $TMPDIR \
--chop-to 100 \
-I {params.block_id_min} \
-R 0 \
-j 0 \
-e 0 \
-l {params.POA_lengths} \
-P {params.POA_params} \
-O 0.001 \
-Y {params.POA_pad_depth} \
-d 0 -D 0 \
-V \
-o {output.gfa}
        '''

rule gffafix:
    input:
        gfa = rules.smoothxg.output['gfa']
    output:
        gfa = 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.gffafix.gfa',
        affixes = 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.gffafix.affixes.tsv.gz'
    threads: 1
    resources:
        mem_mb_per_cpu = 15000,
        runtime = '4h'
    shell:
        '''
        gfaffix {input.gfa} -o {output.gfa} |\
        pigz -p 2 -c > {output.affixes}
        '''

rule vg_path_normalise:
    input:
        fasta = multiext('data/freeze_1/{graph}/{chromosome}.fa.gz','','.fai','.gzi'),
        gfa = rules.gffafix.output['gfa']
    output:
        gfa = 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.vg.gfa'
    params:
        reference = lambda wildcards, input: open(input.fasta[1]).readline().rstrip().split('\t')[0]
    threads: 1
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '24h'
    shell:
        '''
vg convert -t {threads} --packed-out --gfa-in {input.gfa} |\
vg mod -t {threads} --chop 1024 - |\
vg paths -t {threads} -n -Q {params.reference} --xg - |\
vg mod -t {threads} --unchop - |\
vg convert -t {threads} --gfa-out --no-wline - > {output.gfa}
        '''

rule odgi_unchop:
    input:
        gfa = rules.vg_path_normalise.output['gfa'] if config.get('Normalise_path',False) else rules.gffafix.output['gfa']
    output:
        og = 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.unchop.og',
        gfa = 'analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.unchop.gfa'
    threads: 6
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '4h'
    shell:
        '''
odgi build -t {threads} -P -g {input.gfa} -o - -O |\
odgi unchop -P -t {threads} -i - -o - |\
odgi sort -P -p Ygs --temp-dir $TMPDIR -t {threads} -i - -o - |\
tee {output.og} |\
odgi view -i - -g > {output.gfa}
        '''

rule odgi_layout:
    input:
        og = rules.odgi_unchop.output['og']
    output:
        layout = multiext('analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.unchop.lay','','.tsv')
    threads: 6
    resources:
        mem_mb_per_cpu = 8000,
        runtime = '4h'
    shell:
        '''
odgi layout \
-i {input.og} \
 -o {output.layout[0]} \
 -T {output.layout[1]} \
 -t {threads} --temp-dir $TMPDIR -P
        '''

rule odgi_draw:
    input:
        og = rules.odgi_unchop.output['og'],
        layout = rules.odgi_layout.output['layout'][0]
    output:
        image = multiext('analyses/pggb/{graph}/p{p}_s{segment_length}/{chromosome}.k{k}.POA{POA}.unchop.{draw}','.png','.svg')
    params:
        draw_paths = lambda wildcards: '-C -w 20' if wildcards.draw == 'path' else ''
    threads: 2
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '4h'
    shell:
        '''
odgi draw -i {input.og} \
-t {threads} \
--coords-in {input.layout} \
--png {output.image[0]} \
--svg {output.image[1]} \
{params.draw_paths} \
-H 1000
        '''
