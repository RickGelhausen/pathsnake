rule callClusterProfilerRIBO:
    input:
        tsv="functional_analysis/{contrast}/ribo_diffex.tsv",
        rlib="rlib/" + os.path.basename(config["database"])
    output:
        ora=temp(expand("functional_analysis/{{contrast}}/RIBO/ORA/tables/results_GO_{ontology}_{regulation}.tsv", regulation=["up", "down"], ontology=["BP", "MF", "CC"])),
        orasimp=temp(expand("functional_analysis/{{contrast}}/RIBO/ORA/tables/results_GO_simplified_{ontology}_{regulation}.tsv", regulation=["up", "down"], ontology=["BP", "MF", "CC"])),
        gesa=temp(expand("functional_analysis/{{contrast}}/RIBO/GSEA/tables/results_GO_gsea_{ontology}.tsv", ontology=["BP", "MF", "CC"]))
    threads: 1
    conda:
        "../envs/functional.yaml"
    params:
        database=config["database"],
        padjCutoff=config["padjCutoff"],
        log2fcCutoff=config["log2fcCutoff"],
        kegg_org_code=config["keggOrgCode"]
    shell:
        """
        mkdir -p functional_analysis/{wildcards.contrast}/RIBO/;
        pathsnake/scripts/gsea.R -d {params.database} -i {input.tsv} -o functional_analysis/{wildcards.contrast}/RIBO/ -p {params.padjCutoff} -l {params.log2fcCutoff} -k {params.kegg_org_code}
        """

rule callClusterProfilerRNA:
    input:
        tsv="functional_analysis/{contrast}/rna_diffex.tsv",
        rlib="rlib/" + os.path.basename(config["database"])
    output:
        ora=temp(expand("functional_analysis/{{contrast}}/RNA/ORA/tables/results_GO_{ontology}_{regulation}.tsv", regulation=["up", "down"], ontology=["BP", "MF", "CC"])),
        orasimp=temp(expand("functional_analysis/{{contrast}}/RNA/ORA/tables/results_GO_simplified_{ontology}_{regulation}.tsv", regulation=["up", "down"], ontology=["BP", "MF", "CC"])),
        gsea=temp(expand("functional_analysis/{{contrast}}/RNA/GSEA/tables/results_GO_gsea_{ontology}.tsv", ontology=["BP", "MF", "CC"]))
    threads: 1
    conda:
        "../envs/functional.yaml"
    params:
        database=config["database"],
        padjCutoff=config["padjCutoff"],
        log2fcCutoff=config["log2fcCutoff"],
        kegg_org_code=config["keggOrgCode"]
    shell:
        """
        mkdir -p functional_analysis/{wildcards.contrast}/RNA/;
        pathsnake/scripts/gsea.R -d {params.database} -i {input.tsv} -o functional_analysis/{wildcards.contrast}/RNA/ -p {params.padjCutoff} -l {params.log2fcCutoff} -k {params.kegg_org_code}
        """

rule bundleResultsRIBOora:
    input:
        rules.callClusterProfilerRIBO.output.ora
    output:
        "functional_analysis/{contrast}/RIBO/ORA/tables/ORA_GO_results.xlsx"
    conda:
        "../envs/pytools.yaml"
    shell:
        """
        python pathsnake/scripts/bundle_results.py -i {input} -o {output} -m "ORA"
        """

rule bundleResultsRIBOoraSimple:
    input:
        rules.callClusterProfilerRIBO.output.orasimp
    output:
        "functional_analysis/{contrast}/RIBO/ORA/tables/ORA_GO_results_simplified.xlsx"
    conda:
        "../envs/pytools.yaml"
    shell:
        """
        python pathsnake/scripts/bundle_results.py -i {input} -o {output} -m "ORA"
        """

rule bundleResultsRIBOgsea:
    input:
        rules.callClusterProfilerRIBO.output.gesa
    output:
        "functional_analysis/{contrast}/RIBO/GSEA/tables/GSEA_GO_results.xlsx"
    conda:
        "../envs/pytools.yaml"
    shell:
        """
        python pathsnake/scripts/bundle_results.py -i {input} -o {output} -m "GSEA"
        """

rule bundleResultsRNAora:
    input:
        rules.callClusterProfilerRNA.output.ora
    output:
        "functional_analysis/{contrast}/RNA/ORA/tables/ORA_GO_results.xlsx"
    conda:
        "../envs/pytools.yaml"
    shell:
        """
        python pathsnake/scripts/bundle_results.py -i {input} -o {output} -m "ORA"
        """

rule bundleResultsRNAoraSimple:
    input:
        rules.callClusterProfilerRNA.output.orasimp
    output:
        "functional_analysis/{contrast}/RNA/ORA/tables/ORA_GO_results_simplified.xlsx"
    conda:
        "../envs/pytools.yaml"
    shell:
        """
        python pathsnake/scripts/bundle_results.py -i {input} -o {output} -m "ORA"
        """

rule bundleResultsRNAgsea:
    input:
        rules.callClusterProfilerRNA.output.gsea
    output:
        "functional_analysis/{contrast}/RNA/GSEA/tables/GSEA_GO_results.xlsx"
    conda:
        "../envs/pytools.yaml"
    shell:
        """
        python pathsnake/scripts/bundle_results.py -i {input} -o {output} -m "GSEA"
        """