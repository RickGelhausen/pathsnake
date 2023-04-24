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
    make_option(c("-l", "--log2FC_cutoff"), type="double", default=1, help="The log2FC cutoff for the enrichment analysis", metavar="double"),
    make_option(c("-k", "--keggOrgCode"), type="character", default=NULL, help="The KEGG ID of the pathway of interest", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$organism_database) || is.null(opt$input_file) || is.null(opt$output_path) || is.null(opt$keggOrgCode)) {
    print_help(opt_parser)
    stop("Please provide all arguments")
}

organism_database <- opt$organism_database
organism <- strsplit(organism_database, "/")[[1]][2]

library(organism, character.only = TRUE, lib.loc = "rlib")

organism_object <- get(organism)

keytypes(organism_object)

data <- read.table(opt$input_file, header = TRUE, sep = "\t")
geneList <- data$log2FC
names(geneList) <- data$Locus_tag

data_up <- filter(data, log2FC > opt$log2FC_cutoff)# & padj < opt$padj_cutoff)
geneList_up <- data_up$log2FC
names(geneList_up) <- data_up$Locus_tag

data_down <- filter(data, log2FC < -opt$log2FC_cutoff)# & padj < opt$padj_cutoff)
geneList_down <- data_down$log2FC
names(geneList_down) <- data_down$Locus_tag

genes_of_interest_up <- sort(geneList_up, decreasing = TRUE)
genes_of_interest_down <- sort(geneList_down, decreasing = FALSE)

length(geneList)
anyDuplicated(names(geneList))

################################################################################################################################################################
################################################# GO ORA
################################################################################################################################################################
ontologies <- c("BP", "MF", "CC")
regulation <- c("up", "down")
genes_of_interest <- list(genes_of_interest_up, genes_of_interest_down)

dir.create(paste(opt$output_path, "/ORA/tables/", sep=""), recursive=TRUE)
dir.create(paste(opt$output_path, "/ORA/plots/", sep=""), recursive=TRUE)
for (c_ont in ontologies) {
    for (i in 1:length(regulation)) {
        d <- godata(organism, ont=c_ont, keytype="SYMBOL")

        ########## enrichGO ##########
        tryCatch({
            ego <- enrichGO(names(genes_of_interest[[i]]), OrgDb=organism, keyType="SYMBOL", ont = c_ont, pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = FALSE, minGSSize=5, maxGSSize=500)
            results <- data.frame(ego@result)
            write.table(results, file = paste(opt$output_path, "/ORA/tables/results_GO_", c_ont, "_", regulation[i],".tsv", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE)

            #[results$p.adjust <= 0.05, ]
            ego_s <- simplify(ego)
            result_s <- data.frame(ego_s@result)
            write.table(result_s, file = paste(opt$output_path, "/ORA/tables/results_GO_simplified_", c_ont, "_", regulation[i],".tsv", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE)

            # tryCatch({
            #     ego_pt <- pairwise_termsim(ego, semData=d)
            #     #write.table(ego_pt@termsim, file = paste(opt$output_path, "/ORA/tables/pairwise_termsim_", c_ont, "_", regulation[i],".tsv", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE)
            #     #ggsave(paste(opt$output_path, "/ORA/plots/treeplot_ego_", c_ont, "_", regulation[i],".pdf", sep=""), plot=treeplot(ego_pt, fontsize=3), height=16, width=32, units=c("cm"), dpi=600)
            # }, error = function(e) {
            #     print("Not enough terms to plot treeplot:")
            #     print(e)
            # })

            tryCatch({
                ggsave(paste(opt$output_path, "/ORA/plots/dotplot_ORA_GO_", c_ont, "_", regulation[i],".pdf", sep=""), plot=dotplot(ego, showCategory=15), height=20, width=16, units=c("cm"), dpi=600)
            }, error = function(e) {
                print("Not enough terms to plot dotplot:")
                print(e)
            })
        }, error = function(e) {
            df <- data.frame(Doubles=double(),
                 Ints=integer(),
                 Factors=factor(),
                 Logicals=logical(),
                 Characters=character(),
                 stringsAsFactors=FALSE)

            write.table(df, file = paste(opt$output_path, "/ORA/tables/results_GO_", c_ont, "_", regulation[i],".tsv", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE)
            print(e)
        })
    }
}

# ################################################################################################################################################################
# ################################################# KEGG ORA
# ################################################################################################################################################################


dir.create(paste(opt$output_path, "/KEGG/tables/", sep=""), recursive=TRUE)
dir.create(paste(opt$output_path, "/KEGG/plots/", sep=""), recursive=TRUE)

for (i in 1:length(regulation)) {
    kk <- enrichKEGG(gene=names(genes_of_interest[[i]]), organism="mma", pvalueCutoff = 0.05, keyType="kegg")
    results <- data.frame(kk)
    write.table(results, file = paste(opt$output_path, "/KEGG/tables/results_KEGG_", regulation[i],".tsv", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE)
}

################################################################################################################################################################
################################################### GO GSEA
################################################################################################################################################################

hjust_value = 1.0
color_limits = c(0.00, 0.05)  # <-- The cutoff for the color scale for the plot.
count_limits = c(1, 21) # <-- The cutoff for the size of the dots. If < 20 the last dot > 20 dissapears.
break_bin = c(5,10,15,20) # <-- Bins to be shown for the count scale.
radius_range = c(1,11) # <-- Size of the dots.

low_color = "orange" # <-- Color for the dots for the low values
high_color = "blue" # <-- Color for the dots for the high values
na_color = "black" # <-- Color for the dots when no values exist

geneList <- data$log2FC
names(geneList) <- data$Locus_tag
geneList <- sort(geneList, decreasing = TRUE)

dir.create(paste(opt$output_path, "/GSEA/tables/", sep=""), recursive=TRUE)
dir.create(paste(opt$output_path, "/GSEA/plots/", sep=""), recursive=TRUE)
for (c_ont in ontologies) {

    ########## gseGO ##########
    tryCatch({
        gse <- gseGO(geneList=geneList, ont = c_ont, keyType = "SYMBOL", minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = organism, nPermSimple = 100000, eps=0)
        results <- data.frame(gse@result)
        write.table(results, file = paste(opt$output_path, "/GSEA/tables/results_GO_gsea_", c_ont, ".tsv", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE)

        gene_count <- results %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) +1)
        dot_df <- left_join(results, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

        tryCatch({
            df_up <- dot_df[dot_df$enrichmentScore >= 0,]
            df_up_sorted <- df_up %>% arrange(desc(enrichmentScore))
            df_up_sorted <- df_up_sorted[1:10,]

            df_down <- dot_df[dot_df$enrichmentScore < 0,]
            df_down_sorted <- df_down %>% arrange(enrichmentScore)
            df_down_sorted <- df_down_sorted[1:10,]

            print(df_up_sorted)
            print(df_down_sorted)

            p_up <- ggplot(df_up_sorted, aes(x = enrichmentScore, y = fct_reorder(Description, enrichmentScore, .desc = FALSE))) +
                    geom_point(aes(size = count, color = p.adjust)) +
                    scale_colour_gradient(limits=color_limits, low=low_color, high=high_color, na.value=na_color, trans="sqrt") +
                    scale_size_binned(limits=count_limits, breaks=break_bin, range=radius_range) +
                    ylab(NULL) +
                    xlim(0.6, 1.0) +
                    ggtitle("Upregulated") +
                    theme_bw(base_size = 12) +
                    theme(plot.title = element_text(hjust = hjust_value))

            p_down <- ggplot(df_down_sorted, aes(x = enrichmentScore, y = fct_reorder(Description, enrichmentScore, .desc = TRUE))) +
                    geom_point(aes(size = count, color = p.adjust)) +
                    scale_colour_gradient(limits=color_limits, low=low_color, high=high_color, na.value=na_color, trans="sqrt") +
                    scale_size_binned(limits=count_limits, breaks=break_bin, range=radius_range) +
                    ylab(NULL) +
                    scale_x_reverse(limits=c(-0.6,-1.0)) +
                    ggtitle("Downregulated") +
                    theme_bw(base_size = 12) +
                    theme(plot.title = element_text(hjust = hjust_value))

            ggsave(paste(opt$output_path, "/GSEA/plots/dotplot_GSEA_GO_", c_ont, "_up.pdf", sep=""), plot=p_up, height=20, width=16, units=c("cm"), dpi=600)
            ggsave(paste(opt$output_path, "/GSEA/plots/dotplot_GSEA_GO_", c_ont, "_down.pdf", sep=""), plot=p_down, height=20, width=16, units=c("cm"), dpi=600)

        }, error = function(e) {
            print("Not enough terms to plot dotplot:")
            print(e)
        })

        # tryCatch({
        #     gsex <- setReadable(gse, organism, "SYMBOL")
        #     ggsave(paste(opt$output_path, "/GSEA/plots/cnetplot_", c_ont, "_", regulation[i],".pdf", sep=""), plot=cnetplot(gse, foldChange=geneList), height=64, width=64, units=c("cm"), dpi=600)
        # }, error = function(e) {
        #     print("Not enough terms to plot treeplot:")
        #     print(e)
        # })

    }, error = function(e) {
        df <- data.frame(Doubles=double(),
                Ints=integer(),
                Factors=factor(),
                Logicals=logical(),
                Characters=character(),
                stringsAsFactors=FALSE)

        write.table(df, file = paste(opt$output_path, "/GSEA/tables/results_gsea_", c_ont, ".tsv", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE)
        print(e)
    })
}



