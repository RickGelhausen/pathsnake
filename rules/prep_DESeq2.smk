rule prepareEnrichmentInput:
    input:
        xlsx=f"{config['diffexFolder']}/{contrast}_sorted.xlsx"
    output:
        rna="funtional_analysis/{contrast}/rna_diffex.tsv"
    conda:
        "../envs/pytools.yaml"
    params:
        id_column=config["idColumn"],
        type=config["diffexType"]
    shell:
        """
        mkdir -p functional_analysis/{wildcards.contrast}
        pathsnake/prep_diffex_tables.py -x {input.xlsx} -o functional_analysis/{wildcards.contrast} -i {params.id_column} -t {params.type}
        """