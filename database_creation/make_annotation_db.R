library(AnnotationForge)
library("tidyr")


args <- commandArgs(trailingOnly = TRUE)

# Set default values
input_folder <- "."
output_dir <- "databases"
version_num <- "0.1"
maintainer_info <- "Maintainer Name <maintainer@email.com>"
author_info <- "Author Name <author@email.com>"
tax_id_num <- "192952"
genus_name <- "Methanosarcina"
species_name <- "mazei"

# Parse command-line arguments
for (i in seq_along(args)) {
  arg <- args[i]
  next_arg <- args[i + 1]

  if (arg == "--input_folder") {
    input_folder <- next_arg
  } elif (arg == "--output_dir") {
    output_dir <- next_arg
  } elif (arg == "--version") {
    version_num <- next_arg
  } elif (arg == "--maintainer") {
    maintainer_info <- next_arg
  } elif (arg == "--author") {
    author_info <- next_arg
  } elif (arg == "--tax_id") {
    tax_id_num <- next_arg
  } elif (arg == "--genus") {
    genus_name <- next_arg
  } elif (arg == "--species") {
    species_name <- next_arg
  }
}

# Construct file paths
go_file <- file.path(input_folder, "go.csv")
ko_file <- file.path(input_folder, "ko.csv")
refseq_file <- file.path(input_folder, "refseq.csv")
uniprot_file <- file.path(input_folder, "uniprot.csv")

library(AnnotationForge)
library("tidyr")

# Read data files
go_data <- read.table(go_file, sep=",", header=TRUE)
colnames(go_data) <- c("GID", "GO", "EVIDENCE")

ko_data <- read.table(ko_file, sep=",", header=TRUE)
colnames(ko_data) <- c("GID", "KO")

symbol_data <- read.table(refseq_file, sep=",", header=TRUE)
colnames(symbol_data) <- c("GID", "SYMBOL")

uniprot_data <- read.table(uniprot_file, sep=",", header=TRUE)
colnames(uniprot_data) <- c("GID", "UNIPROT")

makeOrgPackage(gene_info=symbol_data, go=go_data, ko=ko_data, uniprot=uniprot_data,
               version=version_num,
               maintainer=maintainer_info,
               author=author_info,
               outputDir=output_dir,
               tax_id=tax_id_num,
               genus=genus_name,
               species=species_name,
               goTable="go")

install.packages("./databases/org.Mmazei.eg.db", repos=NULL)
library("org.Mmazei.eg.db", character.only = TRUE)

keytypes(org.Mmazei.eg.db)