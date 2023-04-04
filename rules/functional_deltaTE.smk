rule callClusterProfilerRIBO:
    input:
        tsv="functional_analysis/{contrast}/ribo_diffex.tsv",
        rlib="rlib/" + os.path.basename(config["database"])
    output:
        bpu="functional_analysis/{contrast}/RIBO/ORA/tables/results_BP_up.tsv",
        bpd="functional_analysis/{contrast}/RIBO/ORA/tables/results_BP_down.tsv",
        mfu="functional_analysis/{contrast}/RIBO/ORA/tables/results_MF_up.tsv",
        mfd="functional_analysis/{contrast}/RIBO/ORA/tables/results_MF_down.tsv",
        ccu="functional_analysis/{contrast}/RIBO/ORA/tables/results_CC_up.tsv",
        ccd="functional_analysis/{contrast}/RIBO/ORA/tables/results_CC_down.tsv",
        gseabp="functional_analysis/{contrast}/RIBO/GSEA/tables/results_gsea_BP.tsv",
        gseamf="functional_analysis/{contrast}/RIBO/GSEA/tables/results_gsea_MF.tsv",
        gseacc="functional_analysis/{contrast}/RIBO/GSEA/tables/results_gsea_CC.tsv"
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
        bpu="functional_analysis/{contrast}/RNA/ORA/tables/results_BP_up.tsv",
        bpd="functional_analysis/{contrast}/RNA/ORA/tables/results_BP_down.tsv",
        mfu="functional_analysis/{contrast}/RNA/ORA/tables/results_MF_up.tsv",
        mfd="functional_analysis/{contrast}/RNA/ORA/tables/results_MF_down.tsv",
        ccu="functional_analysis/{contrast}/RNA/ORA/tables/results_CC_up.tsv",
        ccd="functional_analysis/{contrast}/RNA/ORA/tables/results_CC_down.tsv",
        gseabp="functional_analysis/{contrast}/RNA/GSEA/tables/results_gsea_BP.tsv",
        gseamf="functional_analysis/{contrast}/RNA/GSEA/tables/results_gsea_MF.tsv",
        gseacc="functional_analysis/{contrast}/RNA/GSEA/tables/results_gsea_CC.tsv"
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