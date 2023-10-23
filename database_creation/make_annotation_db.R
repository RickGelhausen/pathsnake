library(AnnotationForge)
library("tidyr")

go_data <- read.table("go.csv", sep=",", header=TRUE)
colnames(go_data) <- c("GID", "GO", "EVIDENCE")

ko_data <- read.table("ko.csv", sep=",", header=TRUE)
colnames(ko_data) <- c("GID", "KO")

symbol_data <- read.table("refseq.csv", sep=",", header=TRUE)
colnames(symbol_data) <- c("GID", "SYMBOL")

uniprot_data <- read.table("uniprot.csv", sep=",", header=TRUE)
colnames(uniprot_data) <- c("GID", "UNIPROT")

dir.create("databases", recursive=TRUE)

makeOrgPackage(gene_info=symbol_data, go=go_data, ko=ko_data, uniprot=uniprot_data,
               version="0.1",
               maintainer="Maintainer Name <maintainer@email.com>",
               author="Maintainer Name <maintainer@email.com>",
               outputDir="databases",
               tax_id="192952",
               genus="Methanosarcina",
               species="mazei",
               goTable="go")

install.packages("./databases/org.Mmazei.eg.db", repos=NULL)
library("org.Mmazei.eg.db", character.only = TRUE)

keytypes(org.Mmazei.eg.db)