# pathsnake

A workflow to perform functional enrichment analysis for prokaryotes using deltaTE or DESeq2 input.
This workflow was developed to be used with the processing pipeline of [HRIBO](https://github.com/RickGelhausen/HRIBO), but will also work with any other correctly formatted input data.

## Creating an annotation database

To run pathsnake which uses the tool [clusterprofiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), an annotation database is needed to perform [Gene Ontology](https://geneontology.org/) Term-based enrichment analysis.

When no curated database is available for an organism (e.g. for [E. coli](https://bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html)), it is possible to create a database based on an existing annotation. In our case, we use the [UniProtKB](https://www.uniprot.org/) to retrieve available GO Term annotation and create an annotation database.

The general idea of this method is to retrieve a table of Gene Identifiers and GO Terms, which can be turned into an SQ-lite Database. This method will work with any table of GO Terms / GeneIDs (and others... It only requires structuring said table correctly).

The scripts and files, we used are available in the `database_creation` subfolder.

:warning: Be advised that the UniprotKB also contains automatic annotations. If you are knowledgeable about ontology in your organism of choice, you can manually improve the input data by pre-filtering it.

### Dependencies

If you want to create your own database, the following dependencies are needed to run the provided scripts:

- python3
- pandas
- r-base
- bioconductor-annotationforge
- bioconductor-go.db

You can install these using conda with the `database.yml`. All tool versions used for our analysis can be found in this file.

```
conda env create -f database.yml

conda activate database_creation
```

### The UniProtKB table

To retrieve the table visit [UniProtKB](https://www.uniprot.org/) and enter an identifier/name to access your organism.

Click on customize columns and make sure all desired columns are visible in the table.
For our example we required: `Entry`, `Gene Names`, `Gene Names (ordered locus)`, `Gene Ontology IDs`, `KEGG`.
Depending on your use-case you might need more or less of these columns.

:warning: Make sure the table only contains the exact organism you want to examine or filter out unwanted identifiers after downloading the table. In our case we filtered for the `METMA' identifier after downloading the table.

### Preparing the UniProt data

After retrieving the table, we need to create a table for each SQL-database table.
To this end we ran `filter_uniprot_table.py`, which extracts the columns from the downloaded table and creates `go.csv`, `ko.csv`, `refseq.csv` and `uniprot.csv`.
These files contain mappings from different identifiers to either GO Terms or KEGG identifiers.

```
python3 filter_uniprot_table.py
```

:warning: There are many ways of doing this step, we just provide one example that we used, but there are many others available. You can use any custom script to filter your data, it is only important to stick to the standards of annotation database Identifiers as can be seen in the example files.

### AnnotationForge

We then use the [AnnotationForge](https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html) R library to create the SQLite database.

To this end adjust the settings in the `make_annotation_db.R` script (i.e. name, email, identifiers for your organism) and run:

```
Rscript make_annotation_db.R
```

This will create an annotation database in the subfolder you specified, with the correct naming convention. This can then be used as an input for `clusterprofiler`.


## Pathsnake

Pathsnake expects as an input an annotation database, as well as an `deltaTE` output folder from `HRIBO`.
It also works with any other correctly formatted differential expression output tables and does not necessarily require `HRIBO`.

### Dependencies

To run the workflow, install [snakemake](https://snakemake.readthedocs.io/en/stable/) all other dependencies will be installed automatically when running `pathsnake`.

You can also install it using conda:

```
conda env create -f snakemake.yml

conda activate snakemake
```

### Preparing the input files

#### database

This file can either be downloaded if a curated version exists for your organism. Otherwise you can create a database yourself using the guide in the previous section.

#### The deltaTE or DESeq output

The differential expression output is required by clusterprofiler to do the analysis.

In fact any differential expression results will work, the difference between both options is that deltaTE provides RNA-seq and RIBO-seq results, while DESeq only provides RNA-seq based results.
As long as they are structured correctly, you can use any differential expression results.

The required format for `deltaTE` is:

- An excel table with the relevant data in the first sheet.
- Naming convention is |treated_condition|-|control|_sorted.xlsx (e.g. B-A_sorted.xlsx)
- It needs to contain the columns: `RIBO_log2FoldChange`, `RIBO_padj`, `RNA_log2FoldChange`, `RNA_padj`

The required format for `deseq2` is:

- An excel table with the relevant data in the first sheet.
- Naming convention is |treated_condition|-|control|_sorted.xlsx (e.g. B-A_sorted.xlsx)
- It needs to contain the columns: `log2FC`, `pvalue_adjusted`

As long as this format is met, any results can be used.

#### Configuration file

Set up the config file:

* Set the path to the differential expression folder.
* Set the idColumn to be used in your provided differential expression tables (e.g. Locus_tag, Old_locus_tag, etc...). :warning: This ID needs to be compatible with the Identifiers used in the annotation database.
* Set the organism code for KEGG.
* Set the organism annotation database for clusterprofiler.

### Run the workflow

Run the workflow by calling:

:warning: Function calls may differ depending on what machine you use and what resources you have available.

```
snakemake -p -k --use-conda -s pathsnake/Snakefile --configfile pathsnake/config.yaml --directory ${PWD} -j 20 --latency-wait 60 &
```

:warning: Snakemake creates a hidden `.snakemake` folder with all the dependencies required to run the tool, make sure to delete it when you no longer need the tool.

## Custom changes

Unfortunately, it is very hard to create universally well-looking plots for every kind input data using `clusterprofiler`.
Known issues might be:

- manual filtering of redundant terms
- very long GO Terms that break the plotting
- custom coloring
- custom setup of conditions
- etc ...

This was also needed for the data used in our publication, therefore we used the script provided in the `custom_scripts` subfolder to customize the appearance of our output results.

We added it for reproducability of our plots and to expose the functions that can be used to customize your plots using `ggplot2`.


## Example

### Database generation

First we will generate the database for methanosarcina mazei.

:exclamation: This step can be skipped, the finished database is in the `example` folder.

1.
