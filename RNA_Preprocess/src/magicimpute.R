#' filename : magicimpute.R
#' Date : 2020-05-17
#' contributor : Yanshuo Chu
#' function: Run magic on seurat obj

suppressMessages({library(optparse)
library(readr)
library(rjson)
library(Seurat)
library(tidyverse)
library(Rmagic)
})

print('---run rmagic---')

option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after normalization",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'seuratobjaftersmooth.rds',
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
## seuobj <- readRDS('../../outs/5GE_VI/original/umap.rds')

###############################################################################
#'                       seuobj after normaliza and scale                    '#
###############################################################################

#' Magic will run PCA itself ##################################################

seuobj.magic <- magic(
    data=seuobj,
    assay = NULL,
    genes = NULL,
    knn = 2, #Note: knn is KA parameter in MAGIC paper
    knn.max = NULL,
    decay = 2, #
    t = 3,
    npca = 200,
    init = NULL,
    t.max = 20,
    knn.dist.method = "euclidean",
    verbose = 1,
    n.jobs = 20
)


saveRDS(seuobj.magic, paste0(opt$out))
