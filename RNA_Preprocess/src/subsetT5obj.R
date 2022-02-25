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
seuobj <- readRDS(opt$data)
##seuobj <- readRDS('../../outs/5GE_VI/original/umap.rds')

#cgoupT <- c(0, 1, 2, 4, 10, 11, 15, 9, 8)
## cgoupT <- c(0,1,2,3,4,9,10,11,12,13,14,17,21)

cgoupT <- c(5)

seuobjT <- subset(seuobj, idents = cgoupT)

saveRDS(seuobjT,  opt$out)
