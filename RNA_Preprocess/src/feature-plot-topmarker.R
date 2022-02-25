#' filename : feature-plot-topmarker.r
#' Date : 2020-04-13
#' contributor : Yanshuo Chu
#'  --

##libraries
suppressMessages({library(optparse)
library(rjson)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
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
colnames(markersTable) <- c("marker", "cluster", "celltype")

clusters <- unique(markersTable$cluster)

for(clt in clusters){

    print("cluster: ")
    print(clt)

    tempMarkerList <- as.vector(markersTable$marker[markersTable$cluster == clt])

    print("marker list: ")
    print(tempMarkerList)
    if (length(tempMarkerList) < 1) {
        next
    }


    gp.list=list()
    for(i in 1: length(tempMarkerList)){
        gp.list = list.append(gp.list, annotate_figure(p = FeaturePlot(seuratobj, reduction = "umap", features = as.vector(tempMarkerList[i]), cols=c("lightgray", "blue", "black")), top = ggpubr::text_grob(label = ct, face="bold", size = 20, color="red")))
    }

    combined.gp <- do.call(ggarrange, c(gp.list, ncol = 2, nrow = 2))
    pdf(paste0(opt$out, "_(cluster: ", clt, ").pdf"), width = param$width, height = param$height)
    print(combined.gp)
    dev.off()
}

