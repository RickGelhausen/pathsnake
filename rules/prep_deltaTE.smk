import os

rule prepareEnrichmentInput:
    input:
        xlsx=os.path.join(config['diffexFolder'], "{contrast}_sorted.xlsx")
    output:
        ribo="functional_analysis/{contrast}/ribo_diffex.tsv",
        rna="functional_analysis/{contrast}/rna_diffex.tsv"
    conda:
        "../envs/pytools.yaml"
    params:
        id_column=config["idColumn"],
        type=config["diffexType"]
    shell:
        """
        mkdir -p functional_analysis/{wildcards.contrast}
        pathsnake/scripts/prep_diffex_tables.py -x {input.xlsx} -o functional_analysis/{wildcards.contrast} -i {params.id_column} -t {params.type}
        """