
#compress output?
rule gbunzip:
    input:
        gbz=rules.cactus_join.output["gbz"]
    output:
        gfa="analyses/cactus/{graph}/graph.gfa"
    shell:
        """
gbunzip --threads {threads} -v -o {output.gfa} {input.gbz}
        """

rule private_graph_nodes:
    input:
        rules.gbunzip.output["gfa"]
    output:
        stats="analyses/cactus/{graph}/graph.private_nodes.stats"
        """
gfa-tools private-nodes -i {input.gfa} -o {output.stats}
        """
