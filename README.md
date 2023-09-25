# pathsnake
A workflow to perform functional enrichment analysis for prokaryotes using deltaTE input.

## Run the workflow

To run the workflow, install snakemake.

Set up the config file:

Set the path to the deltaTE output folder retrieved from HRIBO.
Set the organism code for KEGG.
Set the organism annotation database for clusterprofiler.

Run the workflow by calling:

```
snakemake -p -k --use-conda -s pathsnake/Snakefile --configfile pathsnake/config.yaml --directory ${PWD} -j 20 --latency-wait 60 &
```
