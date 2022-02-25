#' filename : umap-harmony.R
#' Date : 2020-04-23
#' contributor : Yanshuo Chu
#' function: run umap for harmony data

##libraries
suppressMessages({library(optparse)
library(readr)
library(rjson)
library(SeuratData)
library(harmony)
library(Seurat)})
print('---snn clustering---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after normalization",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'snn-harmony.rds',
                help = 'output file name for the r data file [default = %default]',
                metavar = 'character'),
    make_option(c("-c","--param"),
                type = 'character',
                help = 'json file name contain function parameters',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}
if(is.null(opt$param)) {
    print_help(opt_parser)
    stop("json file name (containing parameters) must be provided", call. = F)
}

##load param
param <- fromJSON(file = opt$param)

##Load data
harmony.data <- readRDS(opt$data)

##run snn clustering
harmony.data <- RunUMAP(object = harmony.data,
                     reduction = "harmony",
                     dims = 1:param$npc,
                     min.dist = param$dist)

harmony.data <- FindNeighbors(object = harmony.data, reduction="harmony", dims = 1:param$npc)
snnharmony.obj <- FindClusters(object = harmony.data, resolution = param$res)

saveRDS(snnharmony.obj, file = opt$out)
print('---end---')
