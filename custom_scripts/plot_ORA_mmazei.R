#!/usr/bin/env Rscript

.libPaths( c( .libPaths(), "rlib") )

library(clusterProfiler, lib.loc = "rlib")
library(enrichplot, lib.loc = "rlib")
library(GOSemSim, lib.loc = "rlib")
library(dplyr, lib.loc = "rlib")
library(ggplot2, lib.loc = "rlib")

library(forcats)
library(optparse)
library(stringr)


option_list = list(
    make_option(c("-i", "--input_file"), type="character", default=NULL, help="File containing a column log2FC and an ID/gene symbol", metavar="character"),
    make_option(c("-d", "--organism_database"), type="character", default=NULL, help="Annotation Database for the organism of interest", metavar="character"),
    make_option(c("-o", "--output_path"), type="character", default=NULL, help="The output path containg all output tables and plots.", metavar="character"),
    make_option(c("-p", "--padj_cutoff"), type="double", default=0.05, help="The p-value cutoff for the enrichment analysis", metavar="double"),
    make_option(c("-l", "--log2FC_cutoff"), type="double", default=1, help="The log2FC cutoff for the enrichment analysis", metavar="double")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$organism_database) || is.null(opt$input_file) || is.null(opt$output_path)) {
    print_help(opt_parser)
    stop("Please provide all arguments")
}

organism_database <- opt$organism_database
organism <- strsplit(organism_database, "/")[[1]][2]

library(organism, character.only = TRUE, lib.loc = "rlib")

organism_object <- get(organism)

keytypes(organism_object)

# data <- read.table(opt$input_file, header = TRUE, sep = "\t")
# geneList <- data$log2FC
# names(geneList) <- data$Locus_tag

# data_up <- filter(data, log2FC > opt$log2FC_cutoff & padj < opt$padj_cutoff)
# geneList_up <- data_up$log2FC
# names(geneList_up) <- data_up$Locus_tag

# data_down <- filter(data, log2FC < -opt$log2FC_cutoff & padj < opt$padj_cutoff)
# geneList_down <- data_down$log2FC
# names(geneList_down) <- data_down$Locus_tag

# genes_of_interest_up <- sort(geneList_up, decreasing = TRUE)
# genes_of_interest_down <- sort(geneList_down, decreasing = FALSE)

# length(geneList)
# anyDuplicated(names(geneList))

################################################################################################################################################################
################################################# GO ORA
################################################################################################################################################################
# ontologies <- c("BP")
# regulation <- c("up", "down")
# genes_of_interest <- list(genes_of_interest_up, genes_of_interest_down)

# dir.create(paste(opt$output_path, sep=""), recursive=TRUE)
# for (c_ont in ontologies) {
#     for (i in 1:length(regulation)) {
#         eg <- bitr(genes_of_interest[[i]], fromType = "SYMBOL", toType = "GID", OrgDb = org.Mmazei.eg.db)

#         ########## enrichGO ##########
#         ego <- enrichGO(gene=eg$GID, OrgDb=organism, keyType="GID", ont = c_ont, pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE, minGSSize=5, maxGSSize=500, universe=names(geneList))
#         results <- data.frame(ego@result)
#         write.table(results, file = paste(opt$output_path, "/enrichGO_", c_ont, "_", regulation[i],".tsv", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE)

#         tryCatch({
#             p <- dotplot(ego, showCategory=15)
#             p <- p + scale_color_gradient(limits = c(0, 0.05), low = "red", high = "blue")
#             ggsave(paste(opt$output_path, "/dotplot_ORA_GO_", c_ont, "_", regulation[i],".pdf", sep=""), plot=p, height=20, width=16, units=c("cm"), dpi=600)
#         }, error = function(e) {
#             print("Not enough terms to plot dotplot:")
#             print(e)
#         })
#     }
# }


gene_data <- read.table(opt$input_file, header = TRUE, sep = "\t")
head(gene_data)
nrow(gene_data)
names(gene_data) <- c("Locus_tag", "log2FC", "padj")

geneList <- setNames(gene_data$log2FC, gene_data$Locus_tag)
gene_ids <- names(geneList)
head(gene_ids)

#### UP ####
gene_data_up <- gene_data %>% filter(log2FC > 0 & padj < 0.05)
nrow(gene_data_up)
names(gene_data_up) <- c("Locus_tag", "log2FC", "padj")

geneList_up <- setNames(gene_data_up$log2FC, gene_data_up$Locus_tag)
gene_ids_up <- names(geneList_up)
head(gene_ids_up)

#### DOWN ####
gene_data_down <- gene_data %>% filter(log2FC < 0 & padj < 0.05)
nrow(gene_data_down)
names(gene_data_down) <- c("Locus_tag", "log2FC", "padj")

geneList_down <- setNames(gene_data_down$log2FC, gene_data_down$Locus_tag)
gene_ids_down <- names(geneList_down)
head(gene_ids_down)

eg <- bitr(gene_ids, fromType = "SYMBOL", toType = c("GID", "GO", "UNIPROT", "KO", "EVIDENCE", "GOALL", "ONTOLOGY", "ONTOLOGYALL"), OrgDb = org.Mmazei.eg.db)
head(eg)
# write xlsx file
write.csv(eg, file = "gene_ids.csv")

############ DOWN-REGULATED ############
eg <- bitr(gene_ids_down, fromType = "SYMBOL", toType = "GO", OrgDb = org.Mmazei.eg.db)
write.csv(eg, "downregulated_ribo_gene_symbols.csv")

ego_BP_down <- enrichGO(gene = eg$GO,
                OrgDb = org.Mmazei.eg.db,
                keyType = "GO",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.2,
                ont = "BP",
                readable = TRUE)

# add title to the plot
dotplot(ego_BP_down, title = "Downregulated_BP", showCategory = 15)
# save the plot
ggsave("Downregulated_RIBO_Biological_Processes.png", width = 8, height = 10, units = "in", dpi = 300)

write.csv(ego_BP_down, "Downregulated_RIBO_Biological_Processes.csv")


############ UP-REGULATED ############
eg <- bitr(gene_ids_up, fromType = "SYMBOL", toType = "GO", OrgDb = org.Mmazei.eg.db)
#write csv
write.csv(eg, "upregulated_ribo_gene_symbols.csv")

# Run the enrichment analysis for biological processes
ego_BP_up <- enrichGO(gene = eg$GO,
                OrgDb = org.Mmazei.eg.db,
                keyType = "GO",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.2,
                ont = "BP",
                readable = TRUE)
# add title to the plot
dotplot(ego_BP_up, title = "Upregulated_RIBO_Biological Processes", showCategory = 15)
# save the plot
ggsave("Upregulated_RIBO_Biological_Processes.png", width = 8, height = 10, units = "in", dpi = 300)

write.csv(ego_BP_up, "Upregulated_RIBO_Biological_Processes.csv")
# geneList <- data$log2FC
# names(geneList) <- data$Locus_tag

# data_up <- filter(data, log2FC > opt$log2FC_cutoff & padj < opt$padj_cutoff)
# geneList_up <- data_up$log2FC
# names(geneList_up) <- data_up$Locus_tag

# data_down <- filter(data, log2FC < -opt$log2FC_cutoff & padj < opt$padj_cutoff)
# geneList_down <- data_down$log2FC
# names(geneList_down) <- data_down$Locus_tag

# genes_of_interest_up <- sort(geneList_up, decreasing = TRUE)
# genes_of_interest_down <- sort(geneList_down, decreasing = FALSE)

# length(geneList)
# anyDuplicated(names(geneList))