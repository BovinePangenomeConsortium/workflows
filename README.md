# Bovine Pangenome Consortium workflows

This is the home for pipelines related to the [Bovine Pangenome Consortium](https://bovinepangenomeconsortium.github.io/).
This is a work in progress, but reach out if you have any questions!


### Using the pipeline

Currently, the pipelines are optimised to run on ETH's cluster Euler.
Tools are assumed to be installed, but eventually conda environments will be supported.

We can test the pipeline itself works (without checking the individual jobs work) with

```
git clone git@github.com:BovinePangenomeConsortium/workflows.git
cd .test
snakemake -s ../Snakefile --configfile config.yaml -n
```

### Details

More details on the motivations and details of the pipeline can be found in the [BPC/docs](https://github.com/BovinePangenomeConsortium/docs) repository.
