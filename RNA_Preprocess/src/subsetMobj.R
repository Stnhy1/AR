#' filename : subsetMobj.R
#' Date : 2020-05-30
#' contributor : Yanshuo Chu
#' function: subset monocytes cell


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

cgoupM <- c(17, 18, 7, 14, 19, 21, 6, 12, 5)

seuobjM <- subset(seuobj, idents = cgoupM)

saveRDS(seuobjM,  opt$out)
