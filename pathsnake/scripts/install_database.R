#!/usr/bin/env Rscript

if (!dir.exists("rlib")) {
  dir.create("rlib")
}

.libPaths( c( .libPaths(), "rlib") )

library(optparse)
library(devtools)

devtools::install_github("GuangchuangYu/GOSemSim", lib = "rlib")
devtools::install_github("YuLab-SMU/clusterProfiler", lib = "rlib")

option_list = list(
    make_option(c("-d", "--organism_database"), type="character", default=NULL, help="Annotation Database for the organism of interest", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$organism_database)) {
    print_help(opt_parser)
    stop("Please provide all arguments")
}

organism_database <- opt$organism_database

install.packages(organism_database, repos = NULL, lib="rlib")
