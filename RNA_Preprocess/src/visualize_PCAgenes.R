#' filename : visualize_PCAgenes.R
#' Date : 2020-07-06
#' contributor : Yanshuo Chu
#' function: visualize_PCAgenes

##libraries
suppressMessages({library(optparse)
library(Seurat)
})
print('---visualize embeding---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after runtsne/runumap)",
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

## pdf(file.path(dirname(opt$data), "PCAgenes.pdf"), height = 300,  width = 12)
## print(VizDimLoadings(seuratObj, dims = 1:50, nfeatures = 100,  reduction = "pca"))
## dev.off()

pdf(file.path(dirname(opt$data), "PCAHeatmapGenes.pdf"), height = 180,  width = 12)
print(DimHeatmap(seuratObj, dims = 1:50, cells = 500, nfeatures = 100, balanced = T))
dev.off()

print('---end---')
