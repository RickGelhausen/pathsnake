import os

rule installDatabase:
    output:
        rlib=directory("rlib/" + os.path.basename(config["database"]))
    threads: 1
    conda:
        "../envs/functional.yaml"
    params:
        database=config["database"]
    shell:
        """
        mkdir -p rlib;
        pathsnake/scripts/install_database.R -d {params.database}
        """