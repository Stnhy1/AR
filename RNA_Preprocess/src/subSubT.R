#' filename : subsetobj.R
#' Date : 2020-04-29
#' contributor : Yanshuo Chu
#' function: subset seurat object

##libraries
suppressMessages({library(optparse)
library(readr)
library(rjson)
library(Seurat)
})
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
                default = 'seuratobjwithbatch.rds',
                help = 'output file name for the r data file [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data
seuratObj <- readRDS(opt$data)

Idents(seuratObj) <- seuratObj@meta.data$seurat_clusters

subGroup <- c(16, 19, 20, 21)

subSeuratObj <- subset(seuratObj, idents = subGroup, invert =T)

saveRDS(subSeuratObj,  opt$out)
