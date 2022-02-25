#' filename : addBatchInfo.R
#' Date : 2020-04-17
#' contributor : Yanshuo Chu
#' function: add batch infor to seurat object (metadata)

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
## seuobj <- readRDS('../../outs/5GE_VI/original/umap.rds')

batch1 <- c("A7", "A7B", "A13", "A13B", "A40", "A40B")
batch2 <- c("45", "45B", "A50", "A50B", "A56", "A56B")
batch3 <- c("62", "62B", "63", "63B", "76", "76B")

control <- c("62B", "63B", "65B", "76B", "A7B", "A13B", "A40B", "45B", "A50B", "A56B")


Mono <- c("7", "45", "62", "63")
noAGroup.Mono <- c(Mono, paste0(Mono, "B"))
Group.Mono <- c(noAGroup.Mono, paste0("A", noAGroup.Mono))

OA <- c("13", "50")
noAGroup.OA <- c(OA, paste0(OA, "B"))
Group.OA <- c(noAGroup.OA, paste0("A", noAGroup.OA))

Combo <- c("40", "56", "65", "76")
noAGroup.Combo <- c(Combo, paste0(Combo, "B"))
Group.Combo <- c(noAGroup.Combo, paste0("A", noAGroup.Combo))

md <- seuobj@meta.data
md$batch = "0"
md$control = "SF"
md$GROUP = "NA"

md$batch[md$orig.ident %in% batch1] = "1"
md$batch[md$orig.ident %in% batch2] = "2"
md$batch[md$orig.ident %in% batch3] = "3"

md$GROUP[md$orig.ident %in% Group.Mono] = "Mono"
md$GROUP[md$orig.ident %in% Group.OA] = "OA"
md$GROUP[md$orig.ident %in% Group.Combo] = "Combo"

md$control[md$orig.ident %in% control] = "Blood"

seuobj@meta.data <- md


saveRDS(seuobj, paste0(opt$out))
