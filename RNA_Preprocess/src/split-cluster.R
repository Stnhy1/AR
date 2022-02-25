#'--------------------------------------------------------------
#' filename : split-cluster.R
#' Date : 2021-03-04
#' contributor : Yanshuo Chu
#' function: split-cluster
#'--------------------------------------------------------------

print('<==== split-cluster ====>')

suppressMessages({
    library(optparse)
    library(Seurat)
    library(tidyverse)
    library(ggplot2)
    library(stringr)
})

option_list = list(
    make_option(c("-i","--indata"),
                type = 'character',
                help = 'data.rds',
                metavar = 'character'),
    make_option(c("-o","--outdata"),
                type = 'character',
                help = 'out data',
                metavar = 'character'),
    make_option(c("-c","--cluster"),
                type = 'character',
                help = 'cluster string, seperated by ;',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$indata)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data

seuratObj <- readRDS(opt$indata)
Idents(seuratObj) <- seuratObj$seurat_clusters

clusters <- opt$cluster
clusters <- as.numeric(str_split(clusters, ";")[[1]])

subSeuratObj <- subset(seuratObj, idents = clusters)
saveRDS(subSeuratObj, opt$outdata)


