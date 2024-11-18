# viehgenom

A temporary home for pipelines related to the [Bovine Pangenome Consortium](https://bovinepangenome.github.io/).


## Quality control


We need to get a rough idea of how well assembled these genomes are, but for many samples we only have the assembly itself, not the raw data.
As such, we have to rely on some proxy metrics.
Some should be self-explanatory, but otherwise here is a short description

 + NG50: N50 measure for contiguity assuming a genome size of 3 GB
 + autosome/X/Y single/duplicated/missing copy: How many cetartiodactyla USCOs do we find in single/duplicated/missing states across the different chromosomes
   + we do this to avoid penalising male haplotype-resolved assemblies which **should** be missing many X chromosome-specific USCOs
 + SNPs/InDels: calculated from `minimap2 -c --cs | paftools.js call` using the [ARS-UCD2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002263795.3/) reference genome
 + autosome/X/Y/MT covered: the fraction of each chromosome covered by those `minimap2` alignments
   + we do this as a proxy to identify the sex of the assembly where not easily available and identify obvious misassemblies


## Pangenome processing


## Data analysis

We can also estimate the average gap-compressed identity of the alignments through `cut -f 1,6,13 *.alignment.paf`

