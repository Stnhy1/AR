#' ridge plot
#' 
#' contributor: yanshuo chu

##libraries
suppressMessages({library(optparse)
library(rjson)
library(Seurat)
})
print('---visualize embeding---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after runtsne/runumap)",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'visualization.pdf',
                help = 'output file name for the plot [default = %default]',
                metavar = 'character'),
    make_option(c("-c","--param"),
                type = 'character',
                help = 'json file name contain function parameters',
                metavar = 'character'),
    make_option(c("-m","--marker"),
                type = 'character',
                help = 'marker celltype data frame',
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
seuratobj <- readRDS(opt$data)

markersTable <- read.table(opt$marker, sep = "\t", head = T)
colnames(markersTable)[1:2] <- c("marker", "celltype")
#celltypes <- c("T cell", "B cell", "NK cell")

celltypes <- param$celltypes

pdf(opt$out, width = param$width, height = param$height)
for(ct in celltypes){
    print("cell type: ")
    print(ct)

    tempMarkerList <- unique(as.vector(markersTable$marker[markersTable$celltype == ct]))

    print("marker list: ")
    print(tempMarkerList)

    print(RidgePlot(seuratobj, features = tempMarkerList, ncol=2))
}
dev.off()

