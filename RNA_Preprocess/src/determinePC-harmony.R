#' filename : determinPC-harmony.R
#' Date : 2020-04-23
#' contributor : Yanshuo Chu
#' function: determine harmony pc for umap

print('----determinePC-harmony----')

##libraries
suppressMessages({library(optparse)
library(readr)
library(rjson)
library(Seurat)})

##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after normalization",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'ElbowPlot-harmony.pdf',
                help = 'output file name for the elbow plot [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data
harmony.data <- readRDS(opt$data)

pdf(file.path(opt$out))
VizDimLoadings(object = harmony.data, dims = 1:5, reduction = 'harmony')
DimPlot(object = harmony.data, reduction = 'harmony', group.by = 'orig.ident')
harmony.data <- ProjectDim(object = harmony.data, reduction = 'harmony')
DimHeatmap(object = harmony.data, dims = 1:12, cells = 500, balanced = TRUE,
           reduction = 'harmony')
ElbowPlot(object = harmony.data, ndims = 50, reduction = 'harmony')
dev.off()
print('----end----')
