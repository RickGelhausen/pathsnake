#!/usr/bin/env Rscript

library(optparse)

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
