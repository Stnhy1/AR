#' feature plot
#' 
#' contributor: yanshuo chu

##libraries
suppressMessages({
    library(optparse)
    library(rjson)
    library(rlist)
    library(Seurat)
    library(ggplot2)
    library(RColorBrewer)
    library(ggpubr)
    library(tidyverse)
    ## library(tidymodels)
    ## library(viridis)
    ## library(ggpubr)
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
                metavar = 'character'),
    make_option(c("--control"),
                action="store_true",
                default=FALSE,
                help = 'plot group split by control [default = %default]')
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
seuratObj <- readRDS(opt$data)

markersTable <- read.table(opt$marker, sep = "\t", head = T)
colnames(markersTable)[1:2] <- c("marker", "celltype")


if(any(c("Cd3d", "Cd79a", "Cd79b", "S100a8") %in% rownames(seuratObj))){
    firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        substr(x, 2, nchar(x)) <- tolower(substr(x, 2, nchar(x)))
        x
    }
    markersTable$marker <- firstup(markersTable$marker)
}

#celltypes <- c("T cell", "B cell", "NK cell")

#celltypes <- param$celltypes
celltypes <- unique(markersTable$celltype)

if(opt$control){
    for(ct in celltypes){
        print("cell type: ")
        print(ct)

        tempMarkerList <- as.vector(unique(markersTable$marker[markersTable$celltype == ct]))

        tempMarkerList <- intersect(tempMarkerList, rownames(seuratObj@assays$RNA))

        print("marker list: ")
        print(tempMarkerList)
        if (length(tempMarkerList) < 1) {
            next
        }

        gp.list=list()
        for(i in 1: length(tempMarkerList)){
            gp.list = list.append(gp.list, annotate_figure(p = FeaturePlot(seuratObj, reduction = "umap", features = as.vector(tempMarkerList[i]), cols=c("lightgray", "blue", "black")), top = ggpubr::text_grob(label = ct, face="bold", size = 20, color="red")))
        }

        combined.gp <- do.call(ggarrange, c(gp.list, ncol = 2, nrow = 2))
        pdf(paste0(opt$out, "_(", ct, ").pdf"), width = param$width, height = param$height)
        print(combined.gp)
        dev.off()
    }
}else{
    for(ct in celltypes){
        print("cell type: ")
        print(ct)

        tempMarkerList <- as.vector(unique(markersTable$marker[markersTable$celltype == ct]))

        tempMarkerList <- intersect(tempMarkerList, rownames(seuratObj@assays$RNA))

        print("marker list: ")
        print(tempMarkerList)
        if (length(tempMarkerList) < 1) {
            next
        }

        gp.list=list()
        for(i in 1: length(tempMarkerList)){
            gp.list = list.append(gp.list, annotate_figure(p = FeaturePlot(seuratObj, reduction = "umap", features = as.vector(tempMarkerList[i]), cols=c("lightgray", "blue", "black")), top = ggpubr::text_grob(label = ct, face="bold", size = 20, color="red")))
        }

        combined.gp <- do.call(ggarrange, c(gp.list, ncol = 2, nrow = 2))
        pdf(paste0(opt$out, "_(", ct, ").pdf"), width = param$width, height = param$height)
        print(combined.gp)
        dev.off()
    }
}
